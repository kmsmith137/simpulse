from __future__ import print_function
import zmq
import msgpack
import time
import numpy as np

"""
Python client for the L1 RPC service.

The RPC service is fully asynchronous; each request sends an integer
"token", and all replies to that query include the token.  Thus it is
possible to send multiple requests and then later retrieve the
results.

This client can talk to multiple RPC servers at once.
"""


class AssembledChunk(object):
    """
    This class represents an "assembled chunk" of CHIME/FRB intensity
    data read from a msgpack-format file.

    You probably want to use *read_msgpack_file* to create one of
    these, and then use *decode* to produce arrays of Intensity and Weights.

    Note that rf_pipelines and bonsai assume a frequency ordering of
    high-to-low, so *intensities*[0,0] is the intensity for the first time
    sample and highest frequency (800 MHz).
    """

    def __init__(self, msgpacked_chunk):
        c = msgpacked_chunk
        # print('header', c[0])
        version = c[1]
        assert version in [1, 2]
        if version == 1:
            assert len(c) == 17
        if version == 2:
            assert len(c) == 21
        self.version = version
        # print('version', version)
        compressed = c[2]
        # print('compressed?', compressed)
        compressed_size = c[3]
        # print('compressed size', compressed_size)
        self.beam = c[4]
        self.nupfreq = c[5]
        self.nt_per_packet = c[6]
        self.fpga_counts_per_sample = c[7]
        self.nt_coarse = c[8]
        self.nscales = c[9]
        self.ndata = c[10]
        self.fpga0 = c[11]
        self.fpgaN = c[12]
        self.binning = c[13]

        self.nt = self.nt_coarse * self.nt_per_packet

        scales = c[14]
        offsets = c[15]
        data = c[16]

        # version 2: extra arguments
        self.frame0_nano = 0
        self.nrfifreq = 0
        self.has_rfi_mask = False
        self.rfi_mask = None
        if self.version == 2:
            self.frame0_nano = c[17]
            self.nrfifreq = c[18]
            self.has_rfi_mask = c[19]
            mask = c[20]
            # to numpy
            mask = np.fromstring(mask, dtype=np.uint8)
            mask = mask.reshape((self.nrfifreq, self.nt // 8))
            # Expand mask!
            self.rfi_mask = np.zeros((self.nrfifreq, self.nt), bool)
            for i in range(8):
                self.rfi_mask[:, i::8] = (mask & (1 << i)) > 0

        if compressed:
            import pybitshuffle

            data = pybitshuffle.decompress(data, self.ndata)

        # Convert to numpy arrays
        self.scales = np.fromstring(scales, dtype="<f4")
        self.offsets = np.fromstring(offsets, dtype="<f4")
        self.scales = self.scales.reshape((-1, self.nt_coarse))
        self.offsets = self.offsets.reshape((-1, self.nt_coarse))
        self.data = np.frombuffer(data, dtype=np.uint8)
        self.data = self.data.reshape((-1, self.nt))

    def __str__(self):
        if self.has_rfi_mask:
            h, w = self.rfi_mask.shape
            masked = np.sum(self.rfi_mask == 0)
            rfistr = "yes, %i freqs, %.2f%% masked" % (
                self.nrfifreq,
                (100.0 * masked) / float(h * w),
            )
        else:
            rfistr = "no"
        return "AssembledChunk: beam %i, nt %i, fpga0 %i, rfi %s" % (
            self.beam,
            self.nt,
            self.fpga0,
            rfistr,
        )

    def decode(self):
        # Returns (intensities,weights) as floating-point
        nf = self.data.shape[1]

        # intensities = np.zeros((nf, self.nt), np.float32)
        # weights     = np.zeros((nf, self.nt), np.float32)
        # weights[(self.data > 0) * (self.data < 255)] = 1.

        # print('Data shape:', self.data.shape)
        # print('Scales shape:', self.scales.shape)
        # print('nupfreq:', self.nupfreq)
        # print('nt_per_packet:', self.nt_per_packet)

        intensities = (
            self.offsets.repeat(self.nupfreq, axis=0).repeat(self.nt_per_packet, axis=1)
            + self.data
            * self.scales.repeat(self.nupfreq, axis=0).repeat(
                self.nt_per_packet, axis=1
            )
        ).astype(np.float32)

        weights = ((self.data > 0) * (self.data < 255)) * np.float32(1.0)

        return intensities, weights

    def time_start(self):
        """
        Returns a unix-time (seconds since 1970) value for the start of this
        chunk of data.
        """
        # Nanoseconds per FPGA count
        fpga_nano = 2560
        return 1e-9 * (
            self.frame0_nano + self.fpga_counts_per_sample * fpga_nano * self.fpga0
        )

    def time_end(self):
        """
        Returns a unix-time (seconds since 1970) value for the end of this
        chunk of data.
        """
        # Nanoseconds per FPGA count
        fpga_nano = 2560
        return 1e-9 * (
            self.frame0_nano
            + self.fpga_counts_per_sample * fpga_nano * (self.fpga0 + self.fpgaN)
        )


def read_msgpack_file(fn):
    """ Reads the given *fn* msgpack-formatted CHIME/FRB intensity
    data ("assembled chunk").
    """
    f = open(fn, "rb")
    m = msgpack.unpackb(f.read())
    return AssembledChunk(m)


class WriteChunkReply(object):
    """
    Python version of rpc.hpp : WriteChunks_Reply: the RPC server's
    reply to a WriteChunks request.
    """

    def __init__(self, beam, fpga0, fpgaN, filename, success, err):
        self.beam = beam
        self.fpga0 = fpga0
        self.fpgaN = fpgaN
        self.filename = filename
        self.success = success
        self.error = err
        self.server = None

    def __str__(self):
        s = "WriteChunkReply(beam %i, fpga %i + %i, filename %s" % (
            self.beam,
            self.fpga0,
            self.fpgaN,
            self.filename,
        )
        if self.success:
            s += ", succeeded)"
        else:
            s += ', failed: "%s")' % self.error
        return s

    def __repr__(self):
        return str(self)


class MaskedFrequencies(object):
    """
    The reply to a *get_masked_frequencies* RPC request: a vector of (beam, where) ->
    spectrum mappings giving the fraction of masked samples.
    """

    def __init__(self, msgpack):
        # print('Masked frequencies: got ', msgpack)
        self.histories = {}
        for m in msgpack:
            (beam, where, nt, hist) = m
            print("Where:", where, type(where))
            where = where.decode()
            self.histories[(beam, where)] = np.array(hist).astype(np.float32) / nt


class SummedMaskedFrequencies(object):
    """
    The reply to a *get_summed_masked_frequencies* RPC request: data
    about the fraction of masked samples, by frequency, for a given
    beam, FPGAcounts time range, and place in the RFI chain.
    """

    @staticmethod
    def parse(msgpack):
        # List, one element per beam
        rtn = []
        for m in msgpack:
            mm = SummedMaskedFrequencies(*m)
            rtn.append(mm)
        return rtn

    def __init__(
        self,
        beam,
        fpga_start,
        fpga_end,
        pos_start,
        nt,
        nf,
        nsamples,
        nsamples_unmasked,
        freqs_masked_array,
    ):
        self.beam = beam
        self.fpga_start = fpga_start
        self.fpga_end = fpga_end
        self.pos_start = pos_start
        self.nt = nt
        self.nf = nf
        self.nsamples = nsamples
        self.nsamples_unmasked = nsamples_unmasked
        self.nsamples_masked = nsamples - nsamples_unmasked
        self.freqs_masked = np.array(freqs_masked_array)

    def __str__(self):
        return (
            "SummedMaskedFreq: beam %i, samples %i + %i, nf %i, total masked %i/%i = %.1f %%, frequency histogram: length %i, type %s"
            % (
                self.beam,
                self.pos_start,
                self.nt,
                self.nf,
                self.nsamples_masked,
                self.nsamples,
                100.0 * (self.nsamples_masked) / float(self.nsamples),
                len(self.freqs_masked),
                self.freqs_masked.dtype,
            )
        )


class PacketRate(object):
    """
    The return type for the *get_packet_rate* RPC call.
    """

    def __init__(self, msgpack):
        self.start = msgpack[0]
        self.period = msgpack[1]
        self.packets = msgpack[2]

    def __str__(self):
        return "PacketRate: start %g, period %g, packets: " % (
            self.start,
            self.period,
        ) + str(self.packets)


class InjectData(object):
    """
    Input datatype for the *inject_data* RPC call.
    """

    # matching rf_pipelines_inventory :: inject_data + ch_frb_l1 :: rpc.hpp
    def __init__(self, beam, mode, fpga0, sample_offsets, data):
        """
        *beam*: integer.
        *mode*: 0 = ADD to stream.
        *fpga0*: FPGAcounts of the reference time when these data should be injected.
          The sample times in *sample_offsets* are relative to this FPGAcounts time.
        *sample_offsets*: numpy array of ints, one per frequency.  The intensity sample
          offset (ie, time offset in ~ milliseconds) when data for this frequency should
          be injected.
        *data*: list of numpy arrays, one per frequency; the data to be injected.

        This is a hybrid sparse representation:
        - every frequency is assumed to be represented
        - each frequency can have a different number of samples (including zero),
          and starts at a different time offset.

        Frequency index *f* will get a set of samples injected
        starting at time sample *s0 + sample_offset[f]*, where *s0* is
        *fpga0* converted to a sample number.  The data in *data[f]*
        will be injected.

        Note that this class doesn't care whether data are ordered with the lowest
        frequency first or last in the array.
        """
        # CHIME/FRB values:
        if len(sample_offsets) != 16384:
            print(
                "Warning: expect sample_offsets to have length 16384, got",
                len(sample_offsets),
            )
        if len(data) != 16384:
            print("Warning: expect *data* list to have length 16384, got", len(data))
        # assert(len(sample_offsets) == 16384)
        # assert(len(data) == 16384)
        self.beam = beam
        self.mode = mode
        self.fpga0 = fpga0
        self.sample_offsets = sample_offsets
        self.data = data

    def pack(self, freq_low_to_high):
        # The wire format assumes frequencies are ordered high to low.
        data = self.data
        offsets = self.sample_offsets
        if freq_low_to_high:
            # Flip!!
            data = list(reversed(data))
            offsets = offsets[::-1]

        ndata = np.array([len(d) for d in data]).astype(np.uint16)
        alldata = np.hstack(data).astype(np.float32)
        # msgpack_binary_vector: version 1 packing.
        v = 1
        msg = [
            self.beam,
            self.mode,
            self.fpga0,
            [v, len(offsets), bytes(offsets.astype(np.int32).data)],
            [v, len(ndata), bytes(ndata.data)],
            [v, len(alldata), bytes(alldata.data)],
        ]
        packer = msgpack.Packer(use_single_float=True, use_bin_type=True)
        b = packer.pack(msg)
        return b


class RpcClient(object):
    def __init__(self, servers, context=None, identity=None, debug=False):
        """
        *servers*: a dict of [key, address] entries, where each
          *address* is a string address of an RPC server: eg,
          'tcp://localhost:5555'.

        *context*: if non-None: the ZeroMQ context in which to create
           my socket.

        *identity*: (ignored; here for backward-compatibility)

        NOTE throughout this class that *timeout* is in *milliseconds*!
        """
        if context is None:
            self.context = zmq.Context()
        else:
            self.context = context

        self.do_debug = debug
        self.servers = servers
        self.sockets = {}
        self.rsockets = {}
        for k, v in servers.items():
            self.sockets[k] = self.context.socket(zmq.DEALER)
            self.sockets[k].connect(v)
            self.rsockets[self.sockets[k]] = k

        self.token = 0
        # Buffer of received messages: token->[(message,socket), ...]
        self.received = {}

    def debug(self, *args):
        if not self.do_debug:
            return
        print(*args)

    def get_statistics(self, servers=None, wait=True, timeout=-1):
        """
        Retrieves statistics from each server.  Return value is one
        item per server.  Each item is a list of dictionaries (string
        to int).
        """
        if servers is None:
            servers = self.servers.keys()
        tokens = []
        for k in servers:
            self.token += 1
            req = msgpack.packb(["get_statistics", self.token])
            tokens.append(self.token)
            self.sockets[k].send(req)
        if not wait:
            return tokens
        parts = self.wait_for_tokens(tokens, timeout=timeout)
        # We expect one message part for each token.
        self.debug("get_statistics raw reply:", parts)
        return [msgpack.unpackb(p[0]) if p is not None else None for p in parts]

    def get_packet_rate(
        self, start=None, period=None, servers=None, wait=True, timeout=-1
    ):
        if servers is None:
            servers = self.servers.keys()
        tokens = []

        if start is None:
            start = 0.0
        if period is None:
            period = 0.0

        for k in servers:
            self.token += 1
            req = msgpack.packb(["get_packet_rate", self.token])
            args = msgpack.packb([float(start), float(period)])
            tokens.append(self.token)
            self.sockets[k].send(req + args)
        if not wait:
            return tokens
        parts = self.wait_for_tokens(tokens, timeout=timeout)
        # print('parts:', parts)
        # We expect one message part for each token.
        return [
            PacketRate(msgpack.unpackb(p[0])) if p is not None else None for p in parts
        ]

    def get_packet_rate_history(
        self,
        l0nodes=None,
        start=None,
        end=None,
        period=None,
        servers=None,
        wait=True,
        timeout=-1,
    ):
        if servers is None:
            servers = self.servers.keys()
        tokens = []

        if start is None:
            start = 0.0
        if end is None:
            end = 0.0
        if period is None:
            period = 0.0
        if l0nodes is None:
            l0nodes = ["sum"]

        for k in servers:
            self.token += 1
            req = msgpack.packb(["get_packet_rate_history", self.token])
            args = msgpack.packb([float(start), float(end), float(period), l0nodes])
            tokens.append(self.token)
            self.sockets[k].send(req + args)
        if not wait:
            return tokens
        parts = self.wait_for_tokens(tokens, timeout=timeout)
        # We expect one message part for each token.
        return [msgpack.unpackb(p[0]) if p is not None else None for p in parts]

    def start_fork(
        self, beam_offset, ipaddr, port, beam=0, servers=None, wait=True, timeout=-1
    ):
        """
        beam=0 means send all beams handled by this node.

        If beam != 0, then beam_offset is the destination beam number (not offset)
        """
        if servers is None:
            servers = self.servers.keys()
        tokens = []
        for k in servers:
            self.token += 1
            req = msgpack.packb(["start_fork", self.token])
            args = msgpack.packb([beam, beam_offset, ipaddr, port])
            tokens.append(self.token)
            self.sockets[k].send(req + args)
        if not wait:
            return tokens
        parts = self.wait_for_tokens(tokens, timeout=timeout)
        # We expect one message part for each token.
        return [msgpack.unpackb(p[0]) if p is not None else None for p in parts]

    def stop_fork(
        self, beam_offset, ipaddr, port, beam=0, servers=None, wait=True, timeout=-1
    ):
        if servers is None:
            servers = self.servers.keys()
        tokens = []
        for k in servers:
            self.token += 1
            req = msgpack.packb(["stop_fork", self.token])
            args = msgpack.packb([beam, beam_offset, ipaddr, port])
            tokens.append(self.token)
            self.sockets[k].send(req + args)
        if not wait:
            return tokens
        parts = self.wait_for_tokens(tokens, timeout=timeout)
        # We expect one message part for each token.
        return [msgpack.unpackb(p[0]) if p is not None else None for p in parts]

    def inject_data(self, inj, freq_low_to_high, servers=None, wait=True, timeout=-1):
        """
        Sends data to be injected.  The beam number, time, and data
        are in the *inj* data structure.  Please see the *InjectData*
        class for a detailed explanation of that data structure.

        *freq_low_to_high*: boolean, required.  If True, frequencies in *inj*
        are assumed to be ordered from lowest frequency to highest frequency.
        This is the order returned by simpulse by default.

        (Note that rf_pipelines and bonsai assume high-to-low frequency
        ordering.)

        """
        if servers is None:
            servers = self.servers.keys()
        tokens = []
        for k in servers:
            self.token += 1
            req = msgpack.packb(["inject_data", self.token])
            args = inj.pack(freq_low_to_high)
            print("inject_data argument:", len(args), "bytes")
            tokens.append(self.token)
            self.sockets[k].send(req + args)
        if not wait:
            return tokens
        parts = self.wait_for_tokens(tokens, timeout=timeout)
        # We expect one message part for each token.
        return [msgpack.unpackb(p[0]) if p is not None else None for p in parts]

    def inject_single_pulse(
        self,
        beam,
        sp,
        fpga0,
        nfreq=16384,
        fpga_counts_per_sample=384,
        fpga_period=2.56e-6,
        spectral_modulation=None,
        **kwargs
    ):
        """
        Injects a *simpulse* simulated pulse *sp* into the given *beam*, starting
        at the reference time *fpga0* in FPGAcounts units.
        """
        # NOTE, some approximations here
        t0, t1 = sp.get_endpoints()
        sample_period = fpga_counts_per_sample * fpga_period
        nt = int(np.ceil((t1 - t0) / sample_period))
        t1x = t0 + nt * sample_period
        nsparse = sp.get_n_sparse(t0, t1x, nt)
        print("Pulse time range:", t0, t1x, "NT", nt, "N sparse:", nsparse)
        sparse_data = np.zeros(nsparse, np.float32)
        sparse_i0 = np.zeros(nfreq, np.int32)
        sparse_n = np.zeros(nfreq, np.int32)
        sp.add_to_timestream_sparse(sparse_data, sparse_i0, sparse_n, t0, t1x, nt, 1.0)
        # convert sparse_data into a list of numpy arrays (one per freq)
        data = []
        ntotal = 0
        if spectral_modulation is not None:
            print("Applying spectral modulation to the pulse...")
            if len(spectral_modulation) != nfreq:
                raise RuntimeError(
                    "spectral_modulation, if given, must be an array with length nfreq"
                )
            for i, n in enumerate(sparse_n):
                data.append(sparse_data[ntotal : ntotal + n] * spectral_modulation[i])
                ntotal += n
            # Normalize the fluence to the value from simpulse after modulation
            data_np = np.hstack(data)
            print("Sparse data sum: {}\nPost modulation sum: {}\nsparse_n: {}".format(sparse_data.sum(), data_np.sum(), sparse_n))
            factor = sp.fluence * sp.pulse_nt * sp.nfreq / data_np.sum()
            print("Multiplicative factor to conserve fluence: {}".format(factor))
            for i in range(len(data)):
                for j in range(len(data[i])):
                    data[i][j] *= factor

        else:
            for n in sparse_n:
                data.append(sparse_data[ntotal : ntotal + n])
                ntotal += n
        fpga_offset = int(fpga0 + (t0 / fpga_period))
        injdata = InjectData(beam, 0, fpga_offset, sparse_i0, data)
        # Simpulse orders frequencies low to high.
        freq_low_to_high = True
        return self.inject_data(injdata, freq_low_to_high, **kwargs)

    ## the acq_beams should perhaps be a list of lists of beam ids,
    ## one list per L1 server.
    def stream(
        self,
        acq_name,
        acq_dev="",
        acq_meta="",
        acq_beams=[],
        new_stream=True,
        servers=None,
        wait=True,
        timeout=-1,
    ):
        if servers is None:
            servers = self.servers.keys()
        tokens = []
        for k in servers:
            self.token += 1
            req = msgpack.packb(["stream", self.token])
            # Ensure correct argument types
            acq_name = str(acq_name)
            acq_dev = str(acq_dev)
            acq_meta = str(acq_meta)
            acq_beams = [int(b) for b in acq_beams]
            new_stream = bool(new_stream)
            args = msgpack.packb([acq_name, acq_dev, acq_meta, acq_beams, new_stream])
            tokens.append(self.token)
            self.sockets[k].send(req + args)
        if not wait:
            return tokens
        parts = self.wait_for_tokens(tokens, timeout=timeout)
        # We expect one message part for each token.
        return [msgpack.unpackb(p[0]) if p is not None else None for p in parts]

    def stream_status(self, servers=None, wait=True, timeout=-1):
        if servers is None:
            servers = self.servers.keys()
        tokens = []
        for k in servers:
            self.token += 1
            req = msgpack.packb(["stream_status", self.token])
            tokens.append(self.token)
            self.sockets[k].send(req)
        if not wait:
            return tokens
        parts = self.wait_for_tokens(tokens, timeout=timeout)
        # We expect one message part for each token.
        return [msgpack.unpackb(p[0]) if p is not None else None for p in parts]

    def list_chunks(self, servers=None, wait=True, timeout=-1):
        """
        Retrieves lists of chunks held by each server.
        Return value is one list per server, containing a list of [beam, fpga0, fpga1, bitmask] entries.
        """
        if servers is None:
            servers = self.servers.keys()
        tokens = []
        for k in servers:
            self.token += 1
            req = msgpack.packb(["list_chunks", self.token])
            tokens.append(self.token)
            self.sockets[k].send(req)
        if not wait:
            return tokens
        parts = self.wait_for_tokens(tokens, timeout=timeout)
        # We expect one message part for each token.
        return [msgpack.unpackb(p[0]) if p is not None else None for p in parts]

    def write_chunks(
        self,
        beams,
        min_fpga,
        max_fpga,
        filename_pattern,
        priority=0,
        dm=0.0,
        dm_error=0.0,
        sweep_width=0.0,
        frequency_binning=0,
        need_rfi=None,
        servers=None,
        wait=True,
        timeout=-1,
        waitAll=True,
    ):
        """
        Asks the RPC servers to write a set of chunks to disk.

        *beams*: list of integer beams to write to disk
        *min_fgpa*, *max_fpga*: range of FPGA-counts to write
        *filename_pattern*: printf filename pattern
        *priority*: of writes.

        When requesting a sweep (NOT CURRENTLY IMPLEMENTED!):
        *dm*, *dm_error*: floats, DM and uncertainty of the sweep to request
        *sweep_width*: float, range in seconds to retrieve around the sweep

        *frequency_binning*: int, the factor by which to bin frequency
         data before writing.
        
        *wait*: wait for the initial replies listing the chunks to be written out.
        *waitAll*: wait for servers to reply that all chunks have been written out.
        """
        if servers is None:
            servers = self.servers.keys()
        # Send RPC requests
        tokens = []
        for k in servers:
            self.token += 1
            hdr = msgpack.packb(["write_chunks", self.token])
            request_args = [
                beams,
                min_fpga,
                max_fpga,
                dm,
                dm_error,
                sweep_width,
                frequency_binning,
                filename_pattern,
                priority,
            ]
            if need_rfi in [True, False]:
                request_args.append(need_rfi)
            req = msgpack.packb(request_args)
            tokens.append(self.token)
            self.sockets[k].send(hdr + req)
        if not wait:
            return tokens
        # This will wait for the initial replies from servers, listing
        # the chunks to be written out.
        parts, servers = self.wait_for_tokens(tokens, timeout=timeout, get_sockets=True)
        # We expect one message part for each token.
        chunklists = [msgpack.unpackb(p[0]) if p is not None else None for p in parts]
        # print('Chunklists:', chunklists)
        # print('Servers:', servers)
        if not waitAll:
            # Parse the results into WriteChunkReply objects, and add
            # the .server value.
            results = []
            for chunks, server in zip(chunklists, servers):
                if chunks is None:
                    results.append(None)
                    continue
                rr = []
                for chunk in chunks:
                    res = WriteChunkReply(*chunk)
                    res.server = self.rsockets[server]
                    rr.append(res)
                results.append(rr)
            return results, tokens

        return self.wait_for_all_writes(chunklists, tokens, timeout=timeout)

    def wait_for_all_writes(self, chunklists, tokens, timeout=-1):
        """
        Wait for notification that all writes from a write_chunks
        call have completed.

        *timeout* is in milliseconds.
        """
        results = []
        if timeout > 0:
            t0 = 1000.0 * time.time()
        for token, chunklist in zip(tokens, chunklists):
            if chunklist is None:
                # Did not receive a reply from this server.
                continue
            N = len(chunklist)
            n = 0
            while n < N:
                [p], [s] = self.wait_for_tokens(
                    [token], timeout=timeout, get_sockets=True
                )
                # adjust timeout
                if timeout > 0:
                    tnow = 1000.0 * time.time()
                    timeout = max(0, t0 + timeout - tnow)
                    t0 = tnow
                if p is None:
                    # timed out
                    break
                chunk = msgpack.unpackb(p[0])
                # print('got chunk', chunk)
                n += 1
                # unpack the reply
                [beam, fpga0, fpgaN, filename, success, err] = chunk
                res = WriteChunkReply(beam, fpga0, fpgaN, filename, success, err)
                res.server = self.rsockets[s]
                results.append(res)
        return results

    def get_writechunk_status(self, filename, servers=None, wait=True, timeout=-1):
        """
        Asks the RPC servers for the status of a given filename whose
        write was requested by a write_chunks request.

        Returns (status, error_message).
        """
        if servers is None:
            servers = self.servers.keys()
        # Send RPC requests
        tokens = []
        for k in servers:
            self.token += 1
            hdr = msgpack.packb(["get_writechunk_status", self.token])
            req = msgpack.packb(filename)
            tokens.append(self.token)
            self.sockets[k].send(hdr + req)
        if not wait:
            return tokens
        # Wait for results...
        results = self.wait_for_tokens(tokens, timeout=timeout)
        # We only expect one response per server
        results = [msgpack.unpackb(p[0]) if p is not None else None for p in results]
        return results

    def shutdown(self, servers=None):
        """
        Sends a shutdown request to each RPC server.
        """
        if servers is None:
            servers = self.servers.keys()
        # Send RPC requests
        tokens = []
        for k in servers:
            self.token += 1
            hdr = msgpack.packb(["shutdown", self.token])
            tokens.append(self.token)
            self.sockets[k].send(hdr)

    def start_logging(self, address, servers=None):
        """
        Sends a request to start chlog logging to the given address.
        """
        if servers is None:
            servers = self.servers.keys()
        # Send RPC requests
        tokens = []
        for k in servers:
            self.token += 1
            hdr = msgpack.packb(["start_logging", self.token])
            req = msgpack.packb(address)
            tokens.append(self.token)
            self.sockets[k].send(hdr + req)

    def stop_logging(self, address, servers=None):
        """
        Sends a request to stop chlog logging to the given address.
        """
        if servers is None:
            servers = self.servers.keys()
        # Send RPC requests
        tokens = []
        for k in servers:
            self.token += 1
            hdr = msgpack.packb(["stop_logging", self.token])
            req = msgpack.packb(address)
            tokens.append(self.token)
            self.sockets[k].send(hdr + req)

    def get_masked_frequencies(self, servers=None, wait=True, timeout=-1):
        """
        Deprecated; prefer get_summed_masked_frequencies.

        Returns the most recent mask stats for all beams and places in
        the RFI chain.
        """
        if servers is None:
            servers = self.servers.keys()
        tokens = []
        for k in servers:
            self.token += 1
            req = msgpack.packb(["get_masked_frequencies", self.token])
            # args = msgpack.packb([float(start), float(period)])
            tokens.append(self.token)
            self.sockets[k].send(req)  # + args)
        if not wait:
            return tokens
        parts = self.wait_for_tokens(tokens, timeout=timeout)
        # We expect one message part for each token.
        return [
            MaskedFrequencies(msgpack.unpackb(p[0])) if p is not None else None
            for p in parts
        ]

    def get_summed_masked_frequencies(
        self,
        fpgamin,
        fpgamax,
        beam=-1,
        where="after_rfi",
        servers=None,
        wait=True,
        timeout=-1,
    ):
        """
        Returns summed statistics for the given range of times (in
        FPGAcounts units) from *fpgamin* to *fpgamax*.

        If *beam* is specified, returns only that beam; otherwise, all
        beams.

        *where*: where in the RFI pipeline to retrieve stats from.  In
         production, we typically have only "before_rfi" and
         "after_rfi" available.  This is before and after the L1 RFI
         chain.

        The return value is a list of *SummedMaskedFrequencies* objects.
        """
        fpgamin = int(fpgamin)
        fpgamax = int(fpgamax)
        if servers is None:
            servers = self.servers.keys()
        tokens = []
        for k in servers:
            self.token += 1
            req = msgpack.packb(["get_masked_frequencies_2", self.token])
            args = msgpack.packb([beam, where, fpgamin, fpgamax])
            tokens.append(self.token)
            self.sockets[k].send(req + args)
        if not wait:
            return tokens
        parts = self.wait_for_tokens(tokens, timeout=timeout)
        # We expect one message part for each token.
        return [
            SummedMaskedFrequencies.parse(msgpack.unpackb(p[0]))
            if p is not None
            else None
            for p in parts
        ]

    def get_max_fpga_counts(self, servers=None, wait=True, timeout=-1):
        """
        Returns the largest (ie, last) FPGAcounts values seen at a
        variety of places in the processing pipeline.
        """
        if servers is None:
            servers = self.servers.keys()
        tokens = []
        for k in servers:
            self.token += 1
            req = msgpack.packb(["get_max_fpga_counts", self.token])
            tokens.append(self.token)
            self.sockets[k].send(req)
        if not wait:
            return tokens
        parts = self.wait_for_tokens(tokens, timeout=timeout)
        # We expect one message part for each token.
        return [msgpack.unpackb(p[0]) if p is not None else None for p in parts]

    def _pop_token(self, t, d=None, get_socket=False):
        """
        Pops a message for the given token number *t*, or returns *d*
        if one does not exist.
        """
        try:
            msgs = self.received[t]
        except KeyError:
            return d
        msg = msgs.pop(0)
        # if list of messages for this token is now empty, delete it
        if len(msgs) == 0:
            del self.received[t]
        # self.received actually contains a list of (socket,message) tuples,
        # so drop the socket part if not requested by the caller
        if get_socket:
            return msg
        msg = msg[1]
        return msg

    def _receive(self, timeout=-1):
        """
        Receives >=1 replies from servers.

        *timeout* in milliseconds.  timeout=0 means only read
         already-queued messages.  timeout=-1 means wait forever.

        Returns True if >1 messages were received.
        """

        def _handle_parts(parts, socket):
            hdr = parts[0]
            token = msgpack.unpackb(hdr)
            rest = parts[1:]
            # print('Received token:', token, ':', rest)
            if token in self.received:
                self.received[token].append((socket, rest))
            else:
                self.received[token] = [(socket, rest)]
            # print('Now received parts for token:', self.received[token])

        received = False

        poll = zmq.Poller()
        for k, v in self.sockets.items():
            poll.register(v, zmq.POLLIN)

        # Read all the messages that are waiting (non-blocking read)
        while True:
            events = poll.poll(timeout=0)
            if len(events) == 0:
                break
            for s, e in events:
                # print('Received reply on socket', s)
                parts = s.recv_multipart()
                _handle_parts(parts, s)
                received = True
        if received:
            return True

        if timeout == 0:
            # Don't do a blocking read
            return received

        # Do one blocking read.
        events = poll.poll(timeout)
        for s, e in events:
            # print('Receive()')
            parts = s.recv_multipart()
            _handle_parts(parts, s)
            received = True

        return received

    def wait_for_tokens(self, tokens, timeout=-1, get_sockets=False):
        """
        Retrieves results for the given *tokens*, possibly waiting for
        servers to reply.

        Returns a list of result messages, one for each *token*.

        *timeout* is in milliseconds.
        """

        self.debug("Waiting for tokens (timeout", timeout, "):", tokens)

        results = {}
        if get_sockets:
            sockets = {}
        todo = [token for token in tokens]
        if timeout > 0:
            # note: time.time() is in seconds.
            t0 = 1000.0 * time.time()
        while len(todo):
            self.debug("Still waiting for tokens:", todo)
            done = []
            for token in todo:
                r = self._pop_token(token, get_socket=get_sockets)
                if r is not None:
                    self.debug("Popped token", token, "->", r)
                    done.append(token)
                    if get_sockets:
                        # unpack socket,result tuple
                        socket, r = r
                        sockets[token] = socket
                    results[token] = r
            for token in done:
                todo.remove(token)
            if len(todo):
                self.debug("receive(timeout=", timeout, ")")
                if not self._receive(timeout=timeout):
                    # timed out
                    self.debug("timed out")
                    break
                # adjust timeout
                if timeout > 0:
                    tnow = 1000.0 * time.time()
                    timeout = max(0, t0 + timeout - tnow)
                    t0 = tnow
        self.debug("Results:", results)
        if not get_sockets:
            return [results.get(token, None) for token in tokens]
        return (
            [results.get(token, None) for token in tokens],
            [sockets.get(token, None) for token in tokens],
        )


import threading


class ChLogServer(threading.Thread):
    def __init__(self, addr="127.0.0.1", port=None, context=None):
        super(ChLogServer, self).__init__()
        # I'm a demon thread! Muahahaha
        self.daemon = True

        if context is None:
            self.context = zmq.Context()
        else:
            self.context = context

        self.socket = self.context.socket(zmq.SUB)
        # self.socket.set(zmq.SUBSCRIBE, '')
        self.socket.subscribe("")

        addr = "tcp://" + addr
        if port is None:
            addr += ":*"
            self.socket.bind(addr)
            addr = self.socket.get(zmq.LAST_ENDPOINT)
        else:
            addr += ":" + str(port)
            self.socket.bind(addr)
        print("Bound to", addr)
        self.address = addr

        self.start()

    def run(self):
        while True:
            parts = self.socket.recv_multipart()
            print("Received:", parts)


if __name__ == "__main__":
    import argparse
    import sys

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--shutdown", action="store_true", help="Send shutdown RPC message?"
    )
    parser.add_argument("--log", action="store_true", help="Start up chlog server?")
    parser.add_argument(
        "--fork",
        nargs=4,
        metavar=("<beam>", "<beam offset>", "<dest ip>", "<dest port>"),
        help="Start forking data to the given IP:port with given beam offset.  beam=0 means all beams",
        action="append",
        default=[],
    )
    parser.add_argument(
        "--stop-fork",
        nargs=4,
        metavar=("<beam>", "<beam offset>", "<dest ip>", "<dest port>"),
        help="Stop forking data to the given IP:port with given beam offset.  beam=-1 and destbeam=-1 means stop all forks",
        action="append",
        default=[],
    )
    parser.add_argument(
        "--write",
        "-w",
        nargs=4,
        metavar="x",  # ['<comma-separated beams>', '<minfpga>', '<maxfpga>', '<filename-pattern>'],
        help="Send write_chunks command: <comma-separated beams> <minfpga> <maxfpga> <filename-pattern>",
        action="append",
        default=[],
    )
    parser.add_argument(
        "--need-rfi",
        help="For --write / --awrite, set need_rfi flag.",
        action="store_true",
        default=None,
    )
    parser.add_argument(
        "--awrite",
        nargs=4,
        metavar="x",
        help="Send async write_chunks command: <comma-separated beams> <minfpga> <maxfpga> <filename-pattern>",
        action="append",
        default=[],
    )
    parser.add_argument(
        "--list",
        action="store_true",
        default=False,
        help="Just send list_chunks command and exit.",
    )
    parser.add_argument(
        "--stats",
        action="store_true",
        default=False,
        help="Just request stats and exit.",
    )
    parser.add_argument(
        "--max-fpga",
        action="store_true",
        default=False,
        help="Request max FPGA counts seen at different places in the pipeline.",
    )
    parser.add_argument(
        "--max-fpga-plot",
        type=float,
        help="Request max FPGA counts seen at different places in the pipeline, for N seconds, and make a plot of the results.",
    )
    parser.add_argument("--identity", help="(ignored)")
    parser.add_argument("--stream", help="Stream to files")
    parser.add_argument("--stream-base", help="Stream base directory")
    parser.add_argument("--stream-meta", help="Stream metadata", default="")
    parser.add_argument(
        "--stream-beams",
        action="append",
        default=[],
        help="Stream a subset of beams.  Can be a comma-separated list of integers.  Can be repeated.",
    )

    parser.add_argument(
        "--rate",
        action="store_true",
        default=False,
        help="Send packet rate matrix request",
    )
    parser.add_argument(
        "--rate-history",
        action="store_true",
        default=False,
        help="Send packet rate history request",
    )
    parser.add_argument(
        "--l0",
        action="append",
        default=[],
        help="Request rate history for the list of L0 nodes",
    )
    parser.add_argument(
        "--inject", action="store_true", default=False, help="Inject some data"
    )
    parser.add_argument(
        "--inject-pulse",
        action="store_true",
        default=False,
        help="Inject a simpulse pulse",
    )
    parser.add_argument(
        "--masked-freqs",
        action="store_true",
        default=False,
        help="Send request for masked frequencies history",
    )
    parser.add_argument(
        "--timeout",
        default=10000,
        type=float,
        help="Set timeout (in milliseconds) for all RPCs.  -1 to wait forever.",
    )
    parser.add_argument(
        "ports", nargs="*", help="Addresses or port numbers of RPC servers to contact"
    )
    opt = parser.parse_args()
    args = opt.ports

    if len(args):
        servers = {}
        for i, a in enumerate(args):
            port = None
            try:
                port = int(a)
                port = "tcp://127.0.0.1:%i" % port
            except:
                port = a

            servers[chr(ord("a") + i)] = port
    else:
        servers = dict(a="tcp://127.0.0.1:5555", b="tcp://127.0.0.1:5556")

    print("Sending to servers:", servers)
    client = RpcClient(servers)

    if opt.log:
        logger = ChLogServer()
        addr = logger.address
        client.start_logging(addr)

    doexit = False

    kwa = dict(timeout=opt.timeout)

    if opt.list:
        chunks = client.list_chunks(**kwa)
        print("Received chunk list:", len(chunks))
        for chunklist in chunks:
            for beam, f0, f1, where in chunklist:
                print("  beam %4i, FPGA range %i to %i" % (beam, f0, f1))
        doexit = True

    if opt.stats:
        stats = client.get_statistics(**kwa)
        for s, server in zip(stats, servers.values()):
            print("", server)
            if s is None:
                print("  None")
                continue
            print()
            for d in s:
                keys = d.keys()
                keys.sort()
                for k in keys:
                    print("  ", k, "=", d[k])
                print()
        doexit = True

    if opt.max_fpga:
        fpgas = client.get_max_fpga_counts(**kwa)
        for f, server in zip(fpgas, servers.values()):
            print("", server)
            if f is None:
                print("  None")
                continue
            print(f)
        doexit = True

    if opt.max_fpga_plot:
        import matplotlib

        matplotlib.use("Agg")
        import pylab as plt
        import time
        import fitsio
        from collections import OrderedDict

        t0 = time.time()

        plots = {}

        columns = {}
        times = []

        while True:
            t1 = time.time()
            if t1 - t0 > opt.max_fpga_plot:
                break
            fpgas = client.get_max_fpga_counts(wait=True, **kwa)
            for f, server in zip(fpgas, servers.values()):
                print("", server)
                if f is None:
                    print("  None")
                    continue
                print(f)
                for where, beam, fpga in f:
                    # python3
                    if not isinstance(where, str):
                        where = where.decode()
                    key = (server, beam, where)
                    if not key in plots:
                        plots[key] = []
                    plots[key].append((t1 - t0, fpga))

                    if beam >= 0:
                        key = "beam_%i_%s" % (beam, where)
                    else:
                        key = where
                    key = key.replace(" ", "_")
                    if not key in columns:
                        columns[key] = []
                    # pad missing times
                    npad = len(times) - len(columns[key])
                    if npad:
                        columns[key].extend([0] * npad)
                    columns[key].append(fpga)
            times.append(t1)

            # Sleep to 0.1-second RPC cadence
            t2 = time.time()
            if t2 - t1 < 0.1:
                time.sleep(t1 + 0.1 - t2)

        try:
            fits = fitsio.FITS("max-fpga.fits", "rw", clobber=True)
            fits.write(
                [np.array(c) for c in columns.values()] + [np.array(times)],
                names=list(columns.keys()) + ["time"],
            )
            fits.close()
        except:
            import traceback

            print_exc()

        plt.clf()
        cc = OrderedDict()
        cc["packet_stream"] = ("k", "packets")
        cc["chunk_flushed"] = ("c", "flushed")
        cc["chunk_retrieved"] = ("m", "retrieved")
        cc["before_rfi"] = ("b", "before RFI")
        cc["after_rfi"] = ("g", "after RFI")
        cc["bonsai"] = ("r", "bonsai")
        leg = {}
        for (server, beam, where), plotvals in plots.items():
            print("beam", beam, "where", where)

            label = ""
            # if len(servers) > 1:
            #     label += server + ' '
            # if beam >= 0:
            #     label += 'beam %i ' % beam
            # label += where
            color, label = cc.get(where, (None, where))
            print("color", color, "legend", label)
            dt = np.array([p[0] for p in plotvals])
            fpga = np.array([p[1] for p in plotvals]).astype(float)
            # print('FPGA range', fpga.min(), fpga.max())
            fpga *= 2.56e-6
            I = np.flatnonzero(fpga > 0)
            p = plt.plot(dt[I], fpga[I], ".-", color=color)
            # label=label,
            if not label in leg:
                leg[label] = p[0]

        ll = [k for (c, k) in cc.values() if k in leg]
        lp = [leg[k] for (c, k) in cc.values() if k in leg]
        plt.xlabel("time (s)")
        plt.ylabel("FPGA count (scaled to seconds) (s)")
        plt.legend(lp, ll, loc="upper left")
        plt.savefig("max-fpga.png")

        doexit = True

    if len(opt.fork):
        for beam, beam_offset, ipaddr, port in opt.fork:
            beam = int(beam)
            beam_offset = int(beam_offset)
            port = int(port)
            client.start_fork(beam_offset, ipaddr, port, beam=beam)
        doexit = True

    if len(opt.stop_fork):
        for beam, beam_offset, ipaddr, port in opt.stop_fork:
            beam = int(beam)
            beam_offset = int(beam_offset)
            port = int(port)
            client.stop_fork(beam_offset, ipaddr, port, beam=beam)
        doexit = True

    if opt.variances:
        import matplotlib

        matplotlib.use("Agg")
        import pylab as plt

        results = client.get_bonsai_variances()
        print("Got results: type", type(results))
        print("Got", len(results), "weights and variances")

        freqs = np.linspace(400, 800, 16384)

        plt.clf()
        # one per server
        minvar = 1e3
        for variances in results:
            for b, w, v in variances:
                print(
                    "Beam",
                    b,
                    ":",
                    len(w),
                    "weights, mean",
                    np.mean(w),
                    "and",
                    len(v),
                    "variances, mean",
                    np.mean(v),
                )
                print("First weights:", w[:16])
                print("First variances:", v[:16])
                v = np.array(list(reversed(v)))
                nz, = np.nonzero(v > 0)
                print("v", len(v), "freqs", len(freqs))
                plt.plot(freqs, v, ".", label="Beam %i" % int(b))
                if len(nz):
                    minvar = min(minvar, min(v[nz]))
        plt.yscale("symlog", linthreshy=minvar)
        yl, yh = plt.ylim()
        plt.ylim(-0.1 * minvar, yh)
        plt.xlabel("Frequency (MHz)")
        plt.ylabel("Variances")
        plt.legend()
        plt.title("Bonsai running variance estimates")
        plt.savefig("bonsai-variances.png")

        plt.clf()
        # one per server
        for variances in results:
            for b, w, v in variances:
                w = np.array(list(reversed(w)))
                print("w", len(w), "freqs", len(freqs))
                plt.plot(freqs, w, ".", label="Beam %i" % int(b))
        plt.ylim(-0.1, 1.5)
        plt.ylabel("Weights")
        plt.xlabel("Frequency (MHz)")
        plt.legend()
        plt.title("Bonsai running weight estimates (unmasked fraction)")
        plt.savefig("bonsai-weights.png")

        doexit = True

    if opt.stream:
        beams = []
        for b in opt.stream_beams:
            # Parse possible comma-separate list of strings
            words = b.split(",")
            for w in words:
                beams.append(int(w))
        patterns = client.stream(
            opt.stream, opt.stream_base, opt.stream_meta, beams, **kwa
        )
        print("Streaming to:", patterns)
        doexit = True

    if opt.rate:
        rates = client.get_packet_rate(**kwa)
        print("Received packet rates:")
        for r in rates:
            print("Got rate:", r)
        doexit = True

    if opt.masked_freqs:
        import matplotlib

        matplotlib.use("Agg")
        import pylab as plt

        freqs = client.get_masked_frequencies(**kwa)
        print("Received masked frequencies:")
        for f in freqs:
            print("  ", f)

            for k, v in f.histories.items():
                (beam, where) = k
                hist = v
                print(
                    "Beam",
                    beam,
                    "at",
                    where,
                    ":",
                    hist.shape,
                    hist.dtype,
                    hist.min(),
                    hist.max(),
                )
                plt.clf()
                plt.imshow(
                    hist,
                    interpolation="nearest",
                    origin="lower",
                    vmin=0,
                    vmax=1,
                    cmap="gray",
                    aspect="auto",
                )
                plt.xlabel("Frequency bin")
                plt.ylabel("Time (s)")
                plt.title("Beam %i, %s" % (beam, where))
                plt.savefig("masked-f-%i-%s.png" % (beam, where))
        doexit = True

    if opt.rate_history:
        kwa = {}
        if len(opt.l0):
            kwa.update(l0nodes=opt.l0)
        rates = client.get_packet_rate_history(start=-20, **kwa)
        print("Received packet rate history:")
        print(rates)
        doexit = True

    if len(opt.write):
        for beams, f0, f1, fnpat in opt.write:
            beams = beams.split(",")
            beams = [int(b, 10) for b in beams]
            f0 = int(f0, 10)
            f1 = int(f1, 10)
            kwargs = kwa.copy()
            if opt.need_rfi:
                kwargs.update(need_rfi=True)
            R = client.write_chunks(beams, f0, f1, fnpat, **kwargs)
            print("Results:")
            if R is not None:
                for r in R:
                    print("  ", r)
        doexit = True

    # Async write requests
    if len(opt.awrite):
        for beams, f0, f1, fnpat in opt.awrite:
            beams = beams.split(",")
            beams = [int(b, 10) for b in beams]
            f0 = int(f0, 10)
            f1 = int(f1, 10)
            kwargs = {}
            if opt.need_rfi:
                kwargs.update(need_rfi=True)
            R = client.write_chunks(beams, f0, f1, fnpat, waitAll=False, **kwargs)
            print("Results:")
            if R is not None:
                for r in R:
                    print("  ", r)
        doexit = True

    if opt.inject:
        beam = 10008
        # fpga0 = 52 * 1024 * 384
        # fpga0 = 5 * 1024 * 384
        # fpga0 = 6820000000
        nfreq = 16384
        sample_offsets = np.zeros(nfreq, np.int32)
        data = []
        for f in range(nfreq):
            sample_offsets[f] = int(0.2 * f)
            data.append(100.0 * np.ones(1000, np.float32))
        print(
            "Injecting data spanning",
            np.min(sample_offsets),
            "to",
            np.max(sample_offsets),
            " samples",
        )
        inj = InjectData(beam, 0, fpga0, sample_offsets, data)
        freq_low_to_high = True
        R = client.inject_data(inj, freq_low_to_high, wait=True, **kwa)
        print("Results:", R)

        doexit = True

    if opt.inject_pulse:
        import simpulse

        pulse_nt = 1024
        nfreq = 16384
        freq_lo = 400.0
        freq_hi = 800.0
        # dm = 50.
        dm = 666.0
        sm = 1.0
        width = 0.003
        fluence = 0.1
        spectral_index = -1.0
        undispersed_t = 5.0
        sp = simpulse.single_pulse(
            pulse_nt,
            nfreq,
            freq_lo,
            freq_hi,
            dm,
            sm,
            width,
            fluence,
            spectral_index,
            undispersed_t,
        )
        # fpga0 = 0
        beam = 10008
        fpga0 = 17730 * 1024 * 384
        R = client.inject_single_pulse(beam, sp, fpga0, wait=True, nfreq=nfreq, **kwa)
        print("Results:", R)

        doexit = True

    if doexit:
        sys.exit(0)

    print("get_statistics()...")
    stats = client.get_statistics(**kwa)
    print("Got stats:", stats)

    print("list_chunks()...")
    stats = client.list_chunks(**kwa)
    print("Got chunks:", stats)

    print()
    print("write_chunks()...")

    minfpga = 38600000
    # maxfpga = 38600000
    maxfpga = 48600000

    R = client.write_chunks(
        [77, 78],
        minfpga,
        maxfpga,
        "chunk-beam(BEAM)-chunk(CHUNK)+(NCHUNK).msgpack",
        **kwa
    )
    print("Got:", R)

    for r in R:
        servers = [r.server]
        print("Received", r, "from server:", r.server)
        X = client.get_writechunk_status(r.filename, servers=servers, **kwa)
        X = X[0]
        print("Result:", X)

    print("Bogus result:", client.get_writechunk_status("nonesuch", **kwa))

    chunks, tokens = client.write_chunks(
        [77, 78],
        minfpga,
        maxfpga,
        "chunk2-beam(BEAM)-chunk(CHUNK)+(NCHUNK).msgpack",
        waitAll=False,
        **kwa
    )
    print("Got chunks:", chunks)
    for chlist in chunks:
        if chlist is None:
            continue
        for chunk in chlist:
            print("Chunk:", chunk)
            [R] = client.get_writechunk_status(
                chunk.filename, servers=[chunk.server], **kwa
            )
            print("Status:", R)
    time.sleep(1)
    for chlist in chunks:
        if chlist is None:
            continue
        for chunk in chlist:
            print("Chunk:", chunk)
            [R] = client.get_writechunk_status(
                chunk.filename, servers=[chunk.server], **kwa
            )
            print("Status:", R)

    if opt.log:
        addr = logger.address
        client.stop_logging(addr)

    if opt.shutdown:
        client.shutdown()
        time.sleep(2)
