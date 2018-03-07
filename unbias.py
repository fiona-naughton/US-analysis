import numpy as np
import math

kB_kj = 0.0083144621 # kJ/K/mol
kB_kcal = 0.0019872041 # kCal/K/mol

def calc_probs(method='coord', values=None, Fval=None, x0=None, k=None,
               ener_unit='kJ', T=323):
    k = 1000
    kB = kB_kj if ener_unit is 'kJ' else kB_kcal #!!!!

    if method == 'delta':
        E_spring = 0.5*k*values**2
    elif method == 'force':
        E_spring = 0.5*k*(values/k)**2
    elif method == 'coord':
        E_spring = 0.5*k*(values-x0)**2
    weights = np.exp((-E_spring + Fval)/(kB*T))
    return weights


class UnbiasSystem:
    def __init__(self, Fvals=None, deltas=None, forces=None, coords=None,
                 ks=None, x0s=None, zrange=None, zbins=100):  ##units for Fvals, forces, etc
        if outfilenm == None:
            outfilenm = 'HIST.dat'

        # assume for now I have coords, but get the other two working, as well
        # automatically determine zrange...

        zbin_vals = [zrange[0]+(i+0.5)*(zrange[1]-zrange[0])/zbins for i in range(zbins)]
        self.zbin_vals = zbin_vals
        self.zdata = coords
        win_Evals = [calc_probs(values=np.array(zbin_vals), x0=x0s[wini], Fval=Fvals[wini], k=ks[wini])
                     for wini in range(len(coords))]
        self.win_Evals = win_Evals
        self.zrange = zrange
        self.zbins = zbins
        self.Fvals = Fvals
        self.x0s = x0s
        self.ks = ks

    def unbias(self, data=None, extents=None, numbins=None):
        zbin_vals = self.zbin_vals
        zbins = self.zbins
        coords = self.zdata

        weight_factors = np.zeros((len(zbin_vals),))
        extents = np.array(extents)
        numdatabins = np.array(numbins)

        data = np.array(data)

        hist = np.zeros(tuple([zbins]+numbins))
        weight_factors = np.zeros((zbins))

        for wini, frame_list in enumerate(coords):
            print wini
            win_Evals = self.win_Evals[wini]
            win_count = 0

            weighti = np.zeros((zbins))

            for framei, distval in enumerate(frame_list):
                if distval >= zrange[1] or distval < zrange[0]:
                    continue

                win_count = win_count+1
                dist_bin = int(np.floor(zbins*(distval-zrange[0])/(zrange[1]-zrange[0])))
                data_vals = np.array([d[framei] for d in data[:, wini]])
                data_bins = np.floor(numdatabins*(data_vals-extents[:,0])
                                     /(extents[:,1]-extents[:,0]))

                for bini, binval in enumerate(data_bins):
                    if data_vals[bini] == extents[bini,1]:
                        data_bins[bini] = binval - 1
                        if binval > numdatabins[bini]:
                            data_bins[bini] = binval - numdatabins[bini]

                    data_sel_list = [dist_bin]+[int(k) for k in data_bins]
                    data_sel = tuple(data_sel_list)

                    B_ix = calc_probs(values=distval, x0=self.x0s[wini],
                                     Fval=self.Fvals[wini], k=self.ks[wini])
                    B_ix_cen = calc_probs(values=zbin_vals[dist_bin], x0=self.x0s[wini],
                                     Fval=self.Fvals[wini], k=self.ks[wini])

                    hist[data_sel] = hist[data_sel] + 1
                    weighti[dist_bin] = B_ix_cen  ## fix this

            weight_factors = weight_factors + weighti * win_count

          for disti, dist_level in enumerate(hist):
            if weight_factors[disti] != 0:
              hist[disti] = ((hist[disti]) / weight_factors[disti])
            else:
              hist[disti] = hist[disti]
          hist = np.array(hist)
          self.hist = hist

    def get_hist(self, axes=None, dist_range=None, outfile=None):
        sum_axes = tuple([i for i in range(len(self.hist.shape)) if i not in axes])
        hist = self.hist
        if dist_range is not None:
            dist_bins = [i for i, val in enumerate(self.zbin_vals)
                        if val >= dist_range[0] and val < dist_range[1]]
            hist = hist[dist_bins[0]:dist_bins[-1]+1]
        unbiased_hist = np.sum(hist, axis=sum_axes)
        if outfile is not None:
            np.savetxt(outfile, unbiased_hist)
        return unbiased_hist

    def unbias1D(self, data=None, dist_range=None, outfile=None, err=False):
        zbin_vals = self.zbin_vals
        win_Evals = self.win_Evals
        coords = self.zdata
        zbins = self.zbins

        dist_bin_count = np.zeros((len(zbin_vals,)))
        data = np.array(data)
        len_thingy = [len(data_entry[0]) for indexes, data_entry in np.ndenumerate(data) if isinstance(data_entry, np.ndarray) and data_entry.shape[0] > 0][0]

        dist_bin_weight = np.zeros((len(zbin_vals),len_thingy))
        dist_bin_err_1 = np.zeros((len(zbin_vals),len_thingy))
        dist_bin_err_2 = np.zeros((len(zbin_vals),len_thingy))

        if dist_range is None:
            zrange = self.zrange

        for wini, frame_list in enumerate(coords):
          win_Evals = self.win_Evals[wini]
          win_count = 0
          for framei, distval in enumerate(frame_list):
              if distval >= zrange[1] or distval < zrange[0]:
                  continue

              dist_bin = int(np.floor(zbins*(distval-zrange[0])/(zrange[1]-zrange[0])))
              dist_bin_count[dist_bin] = dist_bin_count[dist_bin] + win_Evals[dist_bin]

              try:
                  dist_bin_weight[dist_bin,:] = dist_bin_weight[dist_bin,:] + win_Evals[dist_bin]*data[wini][framei,:]
                  dist_bin_err_1[dist_bin,:] = dist_bin_err_1[dist_bin,:] + win_Evals[dist_bin]*data[wini][framei,:]*data[wini][framei,:]
                  dist_bin_err_2[dist_bin,:] = dist_bin_err_2[dist_bin,:] + win_Evals[dist_bin]
              except TypeError:
                  continue

        hist = np.zeros((dist_bin_weight.shape[0], dist_bin_weight.shape[1]))

        if err:
            hist = np.zeros((dist_bin_weight.shape[0],dist_bin_weight.shape[1]*2))
        dist_bin_err = np.sqrt((dist_bin_err_1 * dist_bin_err_2 -
                                dist_bin_weight*dist_bin_weight)/
                                (dist_bin_err_2*dist_bin_err_2))

        for j in range(len(zbin_vals)):
            hist[j,:dist_bin_weight.shape[1]] = dist_bin_weight[j,:] / dist_bin_count[j]
        if err:
            hist[:,dist_bin_weight.shape[1]:] = dist_bin_err
        np.savetxt(outfile, hist)
