## Importing modules

# RPY2 is a python package connecting R to python
import rpy2
# Numpy used for data manipulation
import numpy
# Pyplot used for data plotting
import matplotlib



# Secondary imports

# --enable support of Numpy python package
from rpy2.robjects import numpy2ri
# --import R-objects
import rpy2.robjects
# Import R packages to expose all R objects of an R package as Python objects
from rpy2.robjects.packages import importr
# --import plotting tools
import matplotlib.cm as cm



## Class
class CoxFatigueAnalysis:
    
    

    ## Class constructor
    def __init__(self, 
                list_of_covariate_names,        # should be lower-case and no spaces, 
                                                # and not be 'fatigue_life' or 'fatigue_survival'
                list_of_covariate_data_types    # numpy or python data types, e.g. float, int, str
                ):
        
        # Class member variables
        
        # --User defined variables
        self._covariate_names = list_of_covariate_names
        self._covariate_data_type = list_of_covariate_data_types
        self._covariate_number = len(list_of_covariate_names)
        
        # --Basic data structures
        self._fatigue_life = None
        self._fatigue_surival = None
        self._list_of_covariates = []
        self._dataframe_keys = None
        
        # --R libraries as python objects
        self._rbase = None
        self._rsurvival = None
        
        # --R environment objects
        self._robjects = rpy2.robjects
        self._globalenv = self._robjects.globalenv
        
        # --Cox properties
        self._cox_function = None
        self._cox_summary = None
        self._cox_covariates = None
        self._cox_covariates_number = None
        # ----blank lists to hold values 
        self._cox_beta = []
        self._cox_hr = []
        self._cox_hr_min95 = []
        self._cox_hr_max95 = []
        self._cox_se = []
        self._cox_p = []
        
        # --Cox zph properties
        self._cox_zph = None
        # ----blank lists to hold values
        self._cox_zph_rho = []
        self._cox_zph_chisq = []
        self._cox_zph_p = []
        
        
        
        # Expose all R objects of an R package as Python objects
        
        # --import R's "base" package
        self._rbase = importr('base')
        # ----import R's "survival" package
        self._rsurvival = importr('survival')
        # # --import R's "utils" package
        # RUTILS = importr('utils')
    
        
        
        # Setup R-environment (other)
        
        # --automatic conversion from Numpy to R-based data structures
        numpy2ri.activate()
        
        
        
        
    ## Cox Analysis
    def cox_regression(self):
        
        # Creating an R dataframe for survival analysis
    
        # --send data and variables vectors to R environment
        self._globalenv['fatigue_life'] = self._fatigue_life
        self._globalenv['fatigue_survival'] = self._fatigue_survival
        
        # --create dataframe in R and add life and survival
        self._robjects.r('r_dataframe <-' +
                         'data.frame(fatigue_life=fatigue_life, fatigue_survival=fatigue_survival)')
        
        # --add variables to dataframe
        for ii in range(self._covariate_number):
            self._globalenv[self._covariate_names[ii]] = self._list_of_covariates[ii]
            self._robjects.r('r_dataframe${} <- {}'.format(self._covariate_names[ii],self._covariate_names[ii]))
        
        # --create dataframe keys (column names)
        self._dataframe_keys = ['fatigue_life', 'fatigue_survival'] + self._covariate_names 
        
        # --send the list to R environment
        self._globalenv['dataframe_keys'] = self._dataframe_keys
        
        
        
        # Survival and cox functions in R-environment
        
        # --run script that does survival and cox calculations
        # ----script moved to separate file for better code readability
        # ----creates a cox_function object in R environenment
        self._robjects.r.source('r_script_cox.R')
        
        # obtain cox_function from R-environment
        self._cox_function = self._globalenv['cox_function']
        
        
        
        # Recipher cox_function
        
        # --summary of cox_function
        self._cox_summary = self._rbase.summary(self._cox_function)
        cox_coefficients = self._cox_summary[6] # contains coefficient values
        cox_intervals = self._cox_summary[7] # contains HR intervals
        
        # --extract internal covariate names assigned by R
        self._cox_covariates = self._rbase.dimnames(cox_coefficients)[0]
        self._cox_covariates_number = len(self._cox_covariates)
        
        for ii in range(self._cox_covariates_number):
            self._cox_beta.append(cox_coefficients[ii]) #+self._cox_covariates_number*0
            self._cox_hr.append(cox_coefficients[ii+self._cox_covariates_number*1])
            self._cox_se.append(cox_coefficients[ii+self._cox_covariates_number*2])
            self._cox_p.append(cox_coefficients[ii+self._cox_covariates_number*4])
            
            self._cox_hr_min95.append(cox_intervals[ii+self._cox_covariates_number*2])
            self._cox_hr_max95.append(cox_intervals[ii+self._cox_covariates_number*3])
        
    
    
    ## Cox ZPH - test for proportionality assumption
    def cox_zph(self):
        
        # simple check if cox_regression() was run previously.
        if len(self._cox_beta):
            #obtain cox.zph result
            self._cox_zph = self._globalenv['cox_zph']
            
            # recipher cox.zph() function output
            for ii in range(self._cox_covariates_number+1):# +1 for GLOBAL results
                self._cox_zph_rho.append(self._cox_zph[0][ii])
                self._cox_zph_chisq.append(self._cox_zph[0][ii+self._cox_covariates_number+1])
                self._cox_zph_p.append(self._cox_zph[0][ii+(self._cox_covariates_number+1)*2])
        else:
            print('cox_zph(): Please, first run cox_regression()')
    
    
    
    
    ## Survfit calculation function
    def get_cox_survfit_1var(self, 
                    input_covariate, # string of variable to survfit
                    input_covariate_data, # corresponding list of values to survfit
                    constant_covariates, # list of strings of variables to keep constant
                    constant_covariates_values # list of respective constants for each variable
                    ):
                        
        # Sizes of input lists
        input_covariate_number = len(input_covariate_data)
        constant_covariates_number = len(constant_covariates)
        
        
        # Convert input constant_variable to a list of constants
        constant_variable_data = []
        for constant in constant_covariates_values:
            constant_variable_data.append(numpy.array([constant]*input_covariate_number))
        
        
        # Send values to R-environment
        self._globalenv['input_covariate_data'] = numpy.array(input_covariate_data)
        
        
        # Create filter dataframe to use later in survfit
        self._robjects.r('filter_dataframe <- data.frame({}=input_covariate_data)'.format(input_covariate))
        for ii in range(constant_covariates_number):
            self._globalenv['filter'+constant_covariates[ii]] = constant_variable_data[ii]
            self._robjects.r('filter_dataframe${} <- {}'.format(constant_covariates[ii], 'filter'+constant_covariates[ii]))


        # Call survfit function
        self._robjects.r('survfit_function <- survfit(cox_function, newdata=filter_dataframe)')
        survfit_function = self._globalenv['survfit_function']
        
        
        # Extract information from survfit_function
        survfit_time = numpy.array(survfit_function[1])
        survfit_probability = numpy.array(survfit_function[5])
        survfit_survival = numpy.array(survfit_function[4])
        
        return survfit_time, survfit_probability, survfit_survival
    
    
    ## Survfit plot function
    def plot_cox_survival_1var(self,
                            axes, # axes where the data will be plotted
                            input_covariate, # string of variable to survfit
                            input_covariate_data, # corresponding list of values to survfit
                            constant_covariates, # list of strings of variables to keep constant
                            constant_covariates_values, # list of respective constants for each variable
                            axes_format=True # overrides the current formating of axes, default is True
                            ):
                                
        # Sizes of input lists
        input_covariate_number = len(input_covariate_data)
        # constant_covariates_number = len(constant_covariates)
        
        # Call survfit member function
        survfit_time, survfit_probability, survfit_survival = self.get_cox_survfit_1var(input_covariate,
                                                                                input_covariate_data,
                                                                                constant_covariates,
                                                                                constant_covariates_values)
        
        # Add 0,1 point for prettier graphs
        survfit_time = numpy.append(0, survfit_time)
        survfit_survival = numpy.append(0, survfit_survival)
        survfit_probability = numpy.vstack([[1]*survfit_probability.shape[1], survfit_probability])
        
        # Plot all input_covariates
        for ii in range(input_covariate_number):
            axes.step(survfit_time, survfit_probability[:,ii], label=input_covariate_data[ii], where='post')
        
        # Plot censored points
        censored_index = numpy.where(0!=survfit_survival)[0] # 0 indicates failed specimen
        if censored_index.size != 0: # 0 size is returned if all specimens have failed
            censored_time = survfit_time[censored_index]
            censored_probability = survfit_probability[censored_index]
            for ii in range(input_covariate_number):
                axes.plot(censored_time, censored_probability[:,ii], 'x')
        
        # Formatting options
        if axes_format:
            # axes.set_xscale('log')
            axes.set_xlabel('Time')
            axes.set_ylabel('Survival Probability Estimate')
            axes.set_ylim([0,1])
            axes.legend(loc='lower left')




    ## Plot Wohler curve
    def plot_wohler_curve(self, 
                            axes, # axes where the data will be plotted
                            load_covariate_name, # covariate corresponding to load or S on S-N curve
                            load_range, # list of load values [min, max]
                            load_resolution, # set the resolution of the load-axis on wohler curve
                            constant_covariates, # name(s) of covariate(s) to keep constant
                            constant_covariates_values, # respective values
                            axes_format = True, # additional formatting options
                            contour_regions_interval = 0.05, # controls the probablity at which color changes
                            contour_min_threshold = 0.001, # minimum probability threshold for color changes and plotting
                            contour_max_threshold = 0.999, # maximum probability threshold for color changes and plotting
                            colormap = cm.viridis, # selection of matplotlib colormap
                                                   # More info: https://matplotlib.org/examples/color/colormaps_reference.html
                            linecolors = 'k', # color of lines
                            linestyles = ['solid', 'dotted', 'solid', 'dotted', 'solid'], # styles of lines
                            # fill_style = 'contour' #
                            ):
        
        # Create a list corresponding to range and resolution
        load_covariate_data = numpy.arange(load_range[0], load_range[1] + load_resolution, load_resolution)
        
        
        # Run survfit 
        survfit_time, survfit_probability, survfit_survival =   self.get_cox_survfit_1var(
                                                                    load_covariate_name,
                                                                    load_covariate_data,
                                                                    constant_covariates,
                                                                    constant_covariates_values)
        
        
        # Plot returned data
        
        # contour() and contourf() are used to plot the curves. More info at:
        # https://matplotlib.org/api/_as_gen/matplotlib.pyplot.contourf.html
        
        # Control at which probability intervals color changes, default is every 5%
        contourfill_regions = ([0, contour_min_threshold] + # almost 0% probability
                            numpy.arange(contour_regions_interval, 1, contour_regions_interval).tolist() +
                            [contour_max_threshold, 1]) # almost 100% probability
        
        # Fill contours
        axes.contourf(survfit_time, 
                    load_covariate_data, 
                    survfit_probability.transpose(), 
                    contourfill_regions, 
                    cmap=colormap)

        # Draw contour lines at almost 0%, lower and upper 95% bounds, 50% and 100%
        line_contour = axes.contour(survfit_time, 
                    load_covariate_data, 
                    survfit_probability.transpose(), 
                    [contour_min_threshold, 0.05, 0.5, 0.95, contour_max_threshold], 
                    colors=linecolors, 
                    linestyles=linestyles)
        axes.clabel(line_contour, inline=True, fontsize=8, use_clabeltext=True, inline_spacing=1)
        
        
        # Optional axes formating properties
        if axes_format:
            axes.set_xlabel('Time')
            axes.set_ylabel('Load')
        
        # Package raw data for return
        raw_values = [survfit_time, survfit_probability, load_covariate_data]
        
        return [raw_values, line_contour, contourfill_regions]
    
    
    
    
    ## Plot contours (mean, min95, max95, failure, or survival) of 2 curves
    
    #internal function
    @staticmethod
    def __plot_contour(ax, 
                        curve1, 
                        curve2,
                        curve_number, 
                        labels, # labels for plits
                        colors
                        ):
        # Unpack data
        contours = [curve1[1], curve2[1]]
        
        # Extract x,y points
        i=0
        for contour in contours:
            coord = contour.allsegs[curve_number][0]
            x = coord[:,0]
            y = coord[:,1]
            ax.plot(x, y, label=labels[i], color=colors[i])
            i+=1
    
    def compare_curves(self, ax, # pyplot axes
                        curve1, # return of plot_wohler_curve()
                        curve2, # return of plot_wohler_curve()
                        curve_to_plot, # index of the contour to plot
                                    # By default of survival percentage: 
                                                # 'failure' -> contour_min_threshold, 
                                                # 'min95' -> 0.05, 
                                                # 'mean' -> 0.5, 
                                                # 'max95' -> 0.95, 
                                                # 'survival' -> contour_max_threshold 
                        labels=['curve 1', 'curve 2'], # labels for plots
                        curve_colors = ['k','grey'] # colors for curves
                        ):
        cases = {
            'failure' : 0,
            'min95' : 1,
            'mean' : 2,
            'max95' : 3,
            'survival' : 4
        }
        self.__plot_contour(ax, curve1, curve2, cases[curve_to_plot], labels, curve_colors)
    

    ## Printing/console/export functions
    def __repr__(self): #return the string when console print() function is called on class object
        
        # simple check if cox_regression() was run previously.
        if len(self._cox_beta):
            return_string = 'Cox proportional hazards model output:\n'
            
            return_string += 'Cox regression on variables '
            for name in self._cox_covariates:
                return_string += '{}, '.format(name)
            return_string += ':\n\t\t\t beta\t se(beta)\t HR-.95\t HR\t HR+.95\t p\n\n'
            
            # add coefficients
            for ii in range(self._cox_covariates_number):
                return_string += '{}'.format(self._cox_covariates[ii]) + '\t'
                return_string += '{:.5}'.format(self._cox_beta[ii]) + '\t'
                return_string += '{:.5}'.format(self._cox_se[ii]) + '\t'
                return_string += '{:.5}'.format(self._cox_hr_min95[ii]) + '\t'
                return_string += '{:.5}'.format(self._cox_hr[ii]) + '\t'
                return_string += '{:.5}'.format(self._cox_hr_max95[ii]) + '\t'
                return_string += '{:.5}'.format(self._cox_p[ii]) + '\n'
                
            return_string += '\n\n'
    

            # simple check if cox_zph() was run previously.
            if len(self._cox_zph_p):
                return_string += '\nProportionality assumption from cox_zhp() (p>0.05):\n\n'
                return_string += '\t\t\t rho\t chisq\t p\n'
                
                # recipher cox.zph() function output
                for ii in range(self._cox_covariates_number):
                    return_string += '{}:\t {:.4f}\t{:.4f}\t{:.5f}\n'.format(self._cox_covariates[ii], self._cox_zph_rho[ii], self._cox_zph_chisq[ii], self._cox_zph_p[ii])
                return_string += 'GLOBAL:\t {:.4f}\t{:.4f}\t {:.5f}\n'.format(self._cox_zph_rho[self._cox_covariates_number], self._cox_zph_chisq[self._cox_covariates_number], self._cox_zph_p[self._cox_covariates_number])
            
            else:
                return_string += 'No cox_zph() results. Please, run cox_regression()'
                
        else:
            return_string = 'No results. Please, run cox_regression()'
        
        
        return return_string
    
    
    def print_cox_r_output(self):
        
        # simple check if cox_regression() was run previously.
        if len(self._cox_beta):
            return_string = str(self._cox_summary) + '\n\n\n'
            
            # simple check if cox_zph() was run previously.
            if len(self._cox_zph_p):
                return_string += 'Call cox.zph():\n'
                return_string += str(self._cox_zph)
            
            else:
                return_string += 'print_cox_r_output(): No cox_zph() results. Please, run cox_regression()'
        else:
            return_string = 'print_cox_r_output(): No results. Please, run cox_regression()'
        
        print(return_string)



  
    
    ## Importing data from csv
    def import_from_csv(self,
                        filepath, # filepath of the csv file containing data
                        separator=',' # separator used in csv file, default is comma
                        ):
                            
        # import using numpy
        imported_data = numpy.genfromtxt(filepath, delimiter=separator, dtype=object)
        
        # assign member variables
        self._fatigue_life = numpy.array(imported_data[:,0], dtype=float)
        self._fatigue_survival = numpy.array(imported_data[:,1], dtype=int)
        
        # convert covariates to proper data type
        covariates = imported_data[:,2:]
        for ii in range(self._covariate_number):
            self._list_of_covariates.append(numpy.array(covariates[:,ii],
                                                        dtype=self._covariate_data_type[ii]))
         
        
        #TODO check if array of covariates equals to number of names
        #if not prompt that the first x columns are taken
    
    #TODO import_from_excel(
    #TODO import_matrix(
    