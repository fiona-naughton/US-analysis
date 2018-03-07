import numpy as np
import os
from subprocess import check_output


class GrossfieldWHAM:
    wham_units = 'kcal' # the units that wham is compiled in
    wham_command = '~/wham/wham/wham'
    def __init__(self, root='Grossfield-Wham', metafile='metadatafile.dat',
                 timeseriesfiles=None, profilefile='profile.dat',
                 source_data=[], ks=[], x0s=[], temps=[], acs=None,
                 data_type='coord', sim_step=20, sim_start=0, inunits='kcal',
                 stupid=False):
        ## TODO - don't init with lists
        ## assuming type, step and start are the same for all windows!

        if not os.path.exists(root):
            os.makedirs(root)
        if isinstance(timeseriesfiles, list):
            self.timeseries_files = [root+'/'+f for f in timeseriesfiles]
        else:
            if timeseriesfiles is None:
                self.timeseries_files = [root+'/timeseries-{}.dat'.format(i+1)
                                         for i in range(len(source_data))]
            else:
                self.timeseries_files = timesefiesfiles(source_data)
        self.meta_file = root+'/'+metafile
        self.profile_file = root+'/'+profilefile
        self.source_data = source_data
        self.num = len(source_data)
        self.ks = ks
        self.x0s = x0s
        self.temps = temps
        self.sim_step = sim_step
        self.sim_start = sim_start
        self.inunits = inunits
        self.data_type = data_type
        self.acs = acs if acs is not None else [1]*self.num


    def check_min_max(self, minx=None, maxx=None):
         skip_list = []
         max_all = None
         min_all = None
         for i, timeseriesfile in enumerate(self.timeseries_files):
             data = np.genfromtxt(timeseriesfile)
             if len(data) == 0:
                 skip_list.append(i)
                 continue
             mini = data[:,1].min()
             maxi = data[:,1].max()
             if minx is not None:
                 if maxi < minx:
                     skip_list.append(i)
                 else:
                     if maxx is not None:
                         if mini > maxx:
                             skip_list.append(i)
                         else:
                             max_all = maxx
                     else:
                         max_all = (max_all if (max_all is not None and maxi < max_all)
                                   else maxi)
                     min_all = minx
             else:
                minx = (min_all if (min_all is not None and mini > min_all)
                       else mini)
         # TODO - simplify this! + raise a warning when we skip a window

         with open(self.meta_file+'.lim', 'w') as out_file, open(self.meta_file) as in_file:
             for j, line in enumerate(in_file.readlines()):
                 if j not in skip_list:
                     out_file.write(line)
         self.num = self.num - len(skip_list)
         return min_all, max_all


    def run_wham(self, write_input=True, min_x=None, max_x=None, start_time=None,
                 end_time=None, periodicity='',
                 num_bins=200, tol=1e-6, calc_temperature=323,
                 numpad=0, num_MC_trials=200, run_bootstrap=True):
        if write_input:
            self.write_input_files(start_time=start_time, end_time=end_time)
        # periodicity!
        min_x, max_x = self.check_min_max(min_x, max_x)
        run_command = [self.wham_command, periodicity, min_x, max_x, num_bins,
                       tol, calc_temperature, numpad, self.meta_file+'.lim',
                       self.profile_file]
        randSeed = 1 # how to deal with random seed - should make an argument?
        if run_bootstrap:
            run_command = run_command+[num_MC_trials, randSeed]
        check_output(common.list_to_string(run_command), shell=True)
        return self.read_output() # pass energy units here?
        # TODO - optional cleanup of input files


    def write_input_files(self, start_time=None, end_time=None):
        if os.path.exists(self.meta_file):
            os.remove(self.meta_file) ##!!!!
        for i in range(self.num):
            self.add_timeseries_file(start_time=start_time,
                                     end_time=end_time, index=i)

    def add_timeseries_file(self, start_time=None, end_time=None, index=None,
                            data=None, k=None, x0=None, temp=None,
                            timeseries_file=None, correl=None):
        # allow to add individially? could probably remove.
        if index is not None:
            data = self.source_data[index]
            k=self.ks[index]
            x0=self.x0s[index]
            correl=self.acs[index]
            timeseries_file = self.timeseries_files[index]
            #temp=self.temps[i]
        start_time = start_time if start_time is not None else self.sim_start
        start_i = (0 if start_time is None
                  else (start_time-self.sim_start)/self.sim_step)
        end_i = (-1 if end_time is None
                  else (end_time-self.sim_start)/self.sim_step)
        data_slice = data[start_i:end_i]
        # ENERGIES?

        xs = [calc_reaction_coord(line, self.data_type, k, x0)
              for line in data_slice]  ## MAKE SURE DATA IS ITERABLE!!!
        ts = [self.sim_start+self.sim_step*j+start_time for j in range(len(xs))]
        out = np.savetxt(timeseries_file, np.array([ts, xs]).T)
        with open(self.meta_file, 'a') as metafile:
            k_wham_units = convert_energy(k, unitin=self.inunits, unitout=self.wham_units)
            metafile_info = [timeseries_file, x0, k_wham_units, correl]
            ## {energy}
            metafile.write(common.list_to_string(metafile_info)+'\n')
        # warn if all skipped


    def read_output(self, out_units=None):
        profiledata = np.genfromtxt(self.profile_file)

        profiledata[:,1] = convert_energy(profiledata[:,1], unitin=self.wham_units,
                                            unitout=out_units) #temp if kt
        profiledata[:,2] = convert_energy(profiledata[:,2], unitin=self.wham_units,
                                            unitout=out_units) #temp if kt
        with open(self.profile_file) as proffile:
            lines=proffile.readlines()
            try:
                F_vals = [float(l.split()[1]) for l in lines[-self.num:]]
            except ValueError:
                F_vals = [float(l.split()[1]) for l in lines[-self.num+1:]]

        #temp hack cos doesn't know to skip windows on reread; will only work for one here
        F_vals = [convert_energy(F_val, unitin=self.wham_units, unitout=out_units)
                  for F_val in F_vals]
        ## will need to adjust if any dind't make it through...
        return profiledata[:,:3], F_vals





def calc_reaction_coord(value=None, val_type=None, k=None, x0=None):
    def coord(value, k, x0):
        return value
    def delta(value, k, x0):
        if x0 is None:
            raise ValueError('Must provide x0 when using a difference-in-'
                             'reaction-coordinate value')
        return value + x0
    def force(value, k, x0):
        if x0 is None or k is None:
             raise ValueError('Must provide x0 and k when using a force value')
        return -value/k + x0
    options = {'coord': coord, 'delta': delta, 'force': force}
    default = 'coord'
    try:
        calc = options[val_type or default]
    except KeyError:
        raise ValueError("{} is not a valid timeseries data type; valid "
                         "options are: {}".format(val_type, options.keys()))
    if value is not None:
        return calc(value, k, x0)


def convert_energy(value=None, unitin=None, unitout=None, temperature=-1):
    unit_conversions = {'kcal': 1, # kilocalorie, 1 kCal = 1 kCal
                        'kj': 0.239006, # kiloJoule, 1 kJ = 0.239006 kCal
                        'kt': temperature*0.239006 # kT
                        }
    default = 'kcal'
    unitin = unitin if unitin is not None else default
    unitout = unitout if unitout is not None else unitin
    to_kJ = {}
    for unit, name in zip([unitin, unitout], ['in', 'out']):
        try:
            to_kJ[name] = unit_conversions[unit.lower()]
        except KeyError:
            raise ValueError("{} is not a valid energy unit; valid options "
                             "are {}".format(unit, conversions.keys()))
    if value is not None:
        new_val = value * to_kJ['in']/to_kJ['out']
        return new_val
