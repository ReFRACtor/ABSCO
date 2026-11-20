"""Disk-backed storage for the absorption coefficients while a table is being made.

``compute.makeABSCO.calcABSCO`` writes each LBLRTM spectrum into a scratch file as soon as
it is computed, instead of keeping every spectrum in memory.  ``compute.combineVMR`` joins
tables of different mixing ratios without copying them, and ``compute.makeNC`` reads the
result one block of frequencies at a time and writes it to the compressed netCDF file.
Peak memory is then set by the size of one block of frequencies, not by the size of the
table.  The scratch files hold the uncompressed float64 table, so the disk needed in the
run directory is the size of the full table (``absco.ram.table_gb``).

Layout of one scratch file (one band, one mixing-ratio combination): raw little-endian
float64, C order, shape ``(n_layers, n_temperatures, n_frequencies)``, so that one LBLRTM
spectrum is one contiguous row.  Rows that were never written are filled with NaN by
:meth:`BandStore.finish`.
"""

import os

import numpy as np

_DTYPE = np.dtype("<f8")


class BandStore:
    """One band of one table: ``(n_layers, n_temperatures, n_frequencies)`` float64 on disk."""

    def __init__(self, path, n_layers, n_temperatures, n_frequencies):
        self.path = os.fspath(path)
        self.shape = (int(n_layers), int(n_temperatures), int(n_frequencies))
        self.written = np.zeros(self.shape[:2], dtype=bool)
        self._fd = os.open(self.path, os.O_RDWR | os.O_CREAT | os.O_TRUNC, 0o644)
        os.ftruncate(self._fd, int(np.prod(self.shape)) * _DTYPE.itemsize)

    def _offset(self, i_layer, i_temp):
        return (i_layer * self.shape[1] + i_temp) * self.shape[2] * _DTYPE.itemsize

    def write_spectrum(self, i_layer, i_temp, spectrum):
        """Store one spectrum.  A spectrum shorter than the band is padded with NaN; a longer
        one is an error, because the file size was set by the first spectrum."""
        spectrum = np.asarray(spectrum, dtype=_DTYPE)
        n = self.shape[2]
        if spectrum.size > n:
            raise RuntimeError("LBLRTM spectrum has %d points but the band was set up for %d "
                               "by the first spectrum; the runs of one band must not grow"
                               % (spectrum.size, n))
        if spectrum.size < n:
            spectrum = np.concatenate([spectrum, np.full(n - spectrum.size, np.nan, dtype=_DTYPE)])
        os.pwrite(self._fd, np.ascontiguousarray(spectrum).tobytes(), self._offset(i_layer, i_temp))
        self.written[i_layer, i_temp] = True

    def finish(self):
        """Fill the rows no spectrum was written to with NaN."""
        nan_row = np.full(self.shape[2], np.nan, dtype=_DTYPE).tobytes()
        for i_layer, i_temp in zip(*np.nonzero(~self.written)):
            os.pwrite(self._fd, nan_row, self._offset(i_layer, i_temp))
        self.written[:] = True

    def read(self, f0, f1, out):
        """Fill ``out`` (shape ``(n_layers, n_temperatures, f1 - f0)``, C order) with
        frequencies ``f0:f1`` of this band."""
        n_l, n_t, n_f = self.shape
        nbytes = (f1 - f0) * _DTYPE.itemsize
        for i_layer in range(n_l):
            for i_temp in range(n_t):
                row = out[i_layer, i_temp]
                got = os.preadv(self._fd, [memoryview(row).cast("B")],
                                self._offset(i_layer, i_temp) + f0 * _DTYPE.itemsize)
                if got != nbytes:
                    raise RuntimeError("short read from %s" % self.path)

    def remove(self):
        if self._fd is not None:
            os.close(self._fd)
            self._fd = None
        if os.path.exists(self.path):
            os.remove(self.path)


class Table:
    """The absorption coefficients of one mixing-ratio combination, all bands, on disk.

    Reads give arrays of shape ``(frequencies, temperatures, layers)``, as the netCDF
    variable is laid out.  ``bands`` is a list of :class:`BandStore`."""

    def __init__(self, bands):
        self.bands = list(bands)
        self.n_freq_bands = [b.shape[2] for b in self.bands]
        self.edges = np.concatenate([[0], np.cumsum(self.n_freq_bands)]).astype(int)
        n_l, n_t = self.bands[0].shape[:2]
        self.shape = (int(self.edges[-1]), n_t, n_l)

    def read(self, f0, f1):
        n_l, n_t = self.bands[0].shape[:2]
        out = np.empty((f1 - f0, n_t, n_l), dtype=_DTYPE)
        for band, lo, hi in zip(self.bands, self.edges[:-1], self.edges[1:]):
            a, b = max(f0, lo), min(f1, hi)
            if a >= b:
                continue
            buf = np.empty((n_l, n_t, b - a), dtype=_DTYPE)
            band.read(a - lo, b - lo, buf)
            out[a - f0:b - f0] = np.transpose(buf, (2, 1, 0))
        return out

    def remove(self):
        for band in self.bands:
            band.remove()


class StackedTable:
    """Tables of different mixing ratios joined along a new last axis, without copying.

    ``children`` are :class:`Table` or :class:`StackedTable` objects of the same shape."""

    def __init__(self, children):
        self.children = list(children)
        self.shape = tuple(self.children[0].shape) + (len(self.children),)

    def read(self, f0, f1):
        out = np.empty((f1 - f0,) + self.shape[1:], dtype=_DTYPE)
        for k, child in enumerate(self.children):
            out[..., k] = child.read(f0, f1)
        return out

    def remove(self):
        for child in self.children:
            child.remove()
