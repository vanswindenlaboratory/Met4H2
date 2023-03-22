"""Author information: NORCE, Lea Starck, lsta@norceresearch.no 

Disclaimer: NORCE is not in any respect responsible for the use of 
the code or resuls thereof, and assumes no responsibility or 
guarantee for the overall functionality."""

import math
import warnings

import numpy as np
import pandas as pd

from scipy import stats as stats
from scipy.stats import norm

import mcerp as mc

import sympy.physics.units as u
from sympy.physics.units.systems import SI
from sympy import diff


class uncertaintyCalculator:
    """Contains 
        - methods for analytical and numerical computation
        - methods that do monte carlo draws imposing correlations on random input distributions
            - using the Cholesky/Iman-Conover algorithm (code from mcerp).
            - using the slightly modified mcerp code to do SVD when Cholesky fails.
            - mcerp code found at: https://github.com/tisimst/mcerp/blob/master/mcerp/correlate.py under BSD 3-Clause license.
        - method for compiling math expressions.
        - methods for SI unit conversion.
        - methods for collecting computation results in data frames and writing output to excel file.
        - methods that make the calculator easy to run with default settings.
        """

    def __init__(self):
        self.monte_carlo_draws = 10000 
        self.input_dictionary = {}
        self.analysis_results = {}
        self.partial_analysis_results = {}
        self.input_dictionary['computation_method'] =  'lhs_montecarlo'
        self.montecarlo_sensitivity_draws = 50
        
    
    def initialize_calculator(self,input_dictionary=None,SI_conversion=True):
        if input_dictionary is not None:
            self.input_dictionary = input_dictionary
        if not 'list_correlation_matrix' in self.input_dictionary:
            self.set_correlation_matrix(correlation_matrix=[[0,0,0]])

        exist_list = list()
        for elm in self.input_dictionary['list_correlation_matrix']:
            if sum(elm) == 0:
                exist_list.append(False)
            else: exist_list.append(True)
        self.partial_analysis_results['correlations_exist'] = exist_list
        if SI_conversion == True:
            # Convert input units to SI.
            input_keys = ['input_quantities_dictionary','expanded_unc_dictionary']
            self.SI_conversion(input_keys=input_keys)
            # Convert output units to SI.
            output_key = 'out_var_label'
            self.SI_conversion(output_key=output_key) 
        self.calculate_standard_uncertainty()

    def convert_temperature_to_kelvin(self,temp,unit,addition=True):
        if unit == 'C':
            if addition == True: kelvin = temp + 273.15
            else: 
                kelvinfactor = temp
                kelvin = kelvinfactor
        if unit == 'F':
            if addition == True: kelvin = (5/9) * (temp - 32) + 273.15
            else: 
                kelvinfactor = (5/9) * temp
                kelvin = kelvinfactor
        if unit == 'R':
            kelvin = (5/9) * temp
        return kelvin, 'K'

    def get_simple_units(self,unit):
        simple_units = list()
        is_simple_unit = False
        if unit.find("(") == -1:
            unit = unit.replace('(','').replace(')','')
            simple_units.append(unit)
            is_simple_unit = True  
        while is_simple_unit == False:
            simple_unit = unit[unit.find("(")+1:unit.find(")")]
            simple_unit = simple_unit.replace('(','').replace(')','')
            unit = unit[unit.find(")")+1:]
            simple_units.append(simple_unit)
            if unit.find("(") == -1:
                break
        simple_units = [x for x in simple_units if x]

        fahrenheit_in_denominator = False
        if is_simple_unit == False:
            if '°F' in unit:
                f_idx = unit.find('°F')   
                slash_idx = unit.find('/')
                if not slash_idx == -1:
                    first_left_p = unit.find('(',slash_idx)
                    number_p = unit.count('(',first_left_p,f_idx)
                    after_slash = unit[slash_idx:]
                    right_p = slash_idx
                    for _ in range(number_p):
                        idx = after_slash.find(')')
                        right_p = right_p + idx +1
                        after_slash = after_slash[idx+1:]
                    if f_idx > first_left_p: 
                        if f_idx < right_p:
                            fahrenheit_in_denominator = True
        return simple_units, is_simple_unit, fahrenheit_in_denominator

    def convert_simple_units_to_SI(self,input_dict):
        """Converts the input variable units to SI units. 
         Potential pitfalls: 
            - kinds_of_variables needs to be updated so it can take on any conceivable input unit.
            - there may have to be made more ifs to handle derived units or units that are not in the sympy module.
            - The user could mislabel the units, f.ex., use C or c as a unit for temperature, but the 
                code will think C is the unit charge, and c is the speed of light. To avoid confusion refer to the
                Excel file which lists supported units under their assigned categories. Also the 
                unit of the computed output should be compared to the user inputted desired output unit. If they do not 
                match throw an error message. 
            - composite units must be supplied with parentheses according to the pattern given in the example files. 
                Only one fraction is allowed within the composite unit. """ 
        

        # SI unit is always the first element in each list.
        kinds_of_variables = {
            'length':['m','mm','cm','dm','km'],
            'pressure':['Pa','kPa','MPa','bar','bars','psi','atm'],
            'temperature':['K','°F','°C','°R','°Ra'],
            'amount_of_substance':['mol','mols','mole','moles'],
            'time':['s','ns','ms','min','h'],
            'mass':['kg','g','tonn'],
            'el_current':['A'],
            'lum_intensity':['cd'],
            'volume':['m^3','mm^3','cm^3','dm^3','km^3','m**3','mm**3','cm**3','dm**3','km**3','ml','cl','dl','l','mL','cL','dL','L'],
            'area':['m^2','mm^2','cm^2','dm^2','km^2','m**2','mm**2','cm**2','dm**2','km**2'],
            'squared_time':['s^2'],
            'cubed_time':['s^3'],
            'dimension_less':[],
            'no_dim':['','-'],
            'energy':['J','eV']}

        inputs_with_SI_units = {}
        recognized_units = [item for sublist in list(kinds_of_variables.values()) for item in sublist] 
        dim_less = ['-',' ','']
        exception_units = ['°F','°C','°R','°Ra','MPa']

        for label in input_dict:
            inputs_list = list()
            for t in range(len(input_dict[label])):
                unit_list, is_simple_unit, in_denominator = self.get_simple_units(input_dict[label][t][1]) # in denominator refers to whether there is fahrenheit in denominator
                value = input_dict[label][t][0]
                composite_unit_collecter = list()
                composite_value_collecter = list()
                for unit in unit_list:
                    if not unit in recognized_units and not input_dict[label][t][1] in dim_less:
                        print('The input variable \''+label+'\' with units \''+input_dict[label][t][1]+'\' is not supported. Please use input variables with recognized units.') 

                    for kind in kinds_of_variables:
                        # Figure out what sort of variable label is and what its main SI/SI derived unit is.
                        if unit in kinds_of_variables[kind]:
                            label_kind = kind
                            SI_unit = kinds_of_variables[kind][0]
                            # Deal with units that are not in the sympy.physics.units module.
                            if unit == 'min' or unit == 'mins':
                                unit = 'minutes'
                            if unit == 'tonn':
                                unit = 'metric_ton'
                            if unit == '°F':
                                value, unit = self.convert_temperature_to_kelvin(value,'F',addition=is_simple_unit)
                                if is_simple_unit == True: 
                                    if in_denominator == True: value = 1/value
                            if unit == '°C':
                                value, unit = self.convert_temperature_to_kelvin(value,'C',addition=is_simple_unit)
                            if unit == '°R' or unit == '°Ra':
                                value, unit = self.convert_temperature_to_kelvin(value,'R')
                            if unit == 'MPa':
                                value = 1000000*value
                                unit = 'Pa'
                            # If the unit is composite of powers of the same base only then one can get around composite units like so:
                            if label_kind == 'volume' or label_kind == 'area':
                                if '^' in unit or '**' in unit:
                                    power_signs = ['^','**']
                                    fail = 0
                                    for sign in power_signs:
                                        try:
                                            splitted = unit.split(sign)
                                            power = float(splitted[1])
                                            un_powered_unit = splitted[0]
                                            un_powered_SI_unit = SI_unit.split('^')[0]
                                            un_powered_value_with_unit = (input_dict[label][t][0]**(1/power))*getattr(
                                                u, un_powered_unit)
                                            un_powered_new_value_with_unit = u.convert_to(
                                                un_powered_value_with_unit, getattr(u, un_powered_SI_unit))
                                            value = float(un_powered_new_value_with_unit.as_coeff_Mul()[0])**power
                                            unit = (un_powered_new_value_with_unit.as_coeff_Mul()[1])**power
                                        except:
                                            fail = fail + 1
                                            if fail > 1: 
                                                print('Notation was not recognized. Please use \'^\' or \'**\' to represent powers in units, e.g., \'m^3\'.')
                                else:
                                    splitted = SI_unit.split('^')
                                    power = float(splitted[1]) # 3
                                    un_powered_SI_unit = splitted[0] # m
                                    SI_unit = getattr(u, un_powered_SI_unit)**power
                                    value_with_unit = input_dict[label][t][0]*getattr(u, unit)
                                    new_value_with_unit = u.convert_to(value_with_unit, SI_unit)
                                    value = float(new_value_with_unit.as_coeff_Mul()[0])
                                    unit = new_value_with_unit.as_coeff_Mul()[1]
                            else:
                                if not input_dict[label][t][0] == float(0):
                                    if input_dict[label][t][1] in exception_units:
                                        value_with_unit = value*getattr(u,unit)
                                    else:
                                        if not unit in dim_less: 
                                            value_with_unit = input_dict[label][t][0]*getattr(u, unit) 
                                        else:
                                            unit = ''                              
                                else:
                                    value_with_unit = 1*getattr(u, unit)
                                if not unit in dim_less:
                                    new_value_with_unit = u.convert_to(
                                        value_with_unit, getattr(u, SI_unit))
                                if input_dict[label][t][0] == float(0):
                                    value = 0
                                else:
                                    if not unit in dim_less: 
                                        value = float(new_value_with_unit.as_coeff_Mul()[0])
                                if not unit in dim_less:
                                    unit = new_value_with_unit.as_coeff_Mul()[1]
                    if is_simple_unit == False:
                        composite_unit_collecter.append(unit)
                        composite_value_collecter.append(value)
                        value = composite_value_collecter
                        unit = composite_unit_collecter
                    inputs_list.append([value,unit])
            inputs_with_SI_units[label] = inputs_list # This is a nested list so inputs_list[t][0] = [val1,val2], inputs_list[t][0] = [unit1,unit2]
        return inputs_with_SI_units

    def SI_conversion(self,input_keys=None, output_key=None):
        """Potenial pitfall:
            - Display units will need to be updated to handle all possible output units from convert_simple_units_to_SI."""

        output_dict = {}
        display_units = {
            'gas_const':['J/K∙mol','((meter**3.0)*(pascal))/((kelvin)*(mole))'],
            'pressure':['Pa','pascal'],
            'temperature':['K','kelvin'],
            'volume':['m³','meter**3.0'],
            'mass':['kg','kilogram'],
            'density':['kg/m³','(kilogram)/(meter**3.0)'],
            'energy_density':['J/m³','(joule)/(meter**3.0)'],
            'energy_flow':['J/s','(joule)/(second)'],
            'energy':['J','joule'],
            'volume_flow':['m³/s','(meter**3.0)/(second)'],
            'dimension_less':['','-',' '],
            'time':['s','second'],
            'mass_rate':['kg/s','(kilogram)/(second)']
            }

        recognized_units = [item for sublist in list(display_units.values()) for item in sublist]
        
        if input_keys == None: keys = [output_key]
        else: keys = input_keys

        for key in keys:
            if input_keys == None:
                input_dict = {self.input_dictionary['out_var_label']: [[42,self.input_dictionary['out_var_unit']]]}   
            else:
                input_dict = self.input_dictionary[key]
            simple_unit_dict = self.convert_simple_units_to_SI(input_dict)
                
            for var in simple_unit_dict:
                output_list = list()
                for t in range(len(input_dict[var])):
                
                    result = input_dict[var][t][0]
                    if isinstance(simple_unit_dict[var][t][0],list):
                        composite_unit = input_dict[var][t][1]
                        old_units, _ , _ = self.get_simple_units(composite_unit)
                        part_no = 0
                        for part in simple_unit_dict[var][t][0]:
                            factor = part/input_dict[var][t][0]
                            if composite_unit.find('/') < composite_unit.find(old_units[part_no]):
                                result = result/factor
                            else:
                                result = factor*result
                            composite_unit = composite_unit.replace(
                                old_units[part_no],str(simple_unit_dict[var][t][1][part_no]))
                            part_no = part_no + 1
                        
                    else: composite_unit = simple_unit_dict[var][t][1]
                    for display_unit in display_units: 
                        if not str(composite_unit) in recognized_units:
                            print('The input variable \''+var+'\' with units \''+str(composite_unit)+'\' is not supported.')
                        else:
                            if str(composite_unit) in display_units[display_unit]:
                                if isinstance(simple_unit_dict[var][t][0],list):
                                        output_list.append([result,display_units[display_unit][0]])
                                else: output_list.append(
                                    [simple_unit_dict[var][t][0],display_units[display_unit][0]])
                output_dict[var] = output_list.copy()
            if key != 'out_var_label':
                self.input_dictionary[key] = output_dict.copy()
            else: self.input_dictionary['out_var_unit'] = output_dict[self.input_dictionary['out_var_label']][0][1]
        
    def correlate_SVD(self,params,correlation_matrix):
        # This method was found at https://github.com/tisimst/mcerp/blob/master/mcerp/correlate.py
        # only slightly modified to accomodate SVD when Cholesky decomposition fails with negative eigenvalues.
        """
        Force a correlation matrix on a set of statistically distributed objects.
        This function works on objects in-place.
    
        Parameters
        ----------
        params : array
            An array of of uv objects.
        correlation_matrix : 2d-array
            The correlation matrix to be imposed
        """
        # Make sure all inputs are compatible
        assert all(
            [isinstance(param, mc.UncertainFunction) for param in params]
        ), 'All inputs to "correlate" must be of type "UncertainFunction"'

        # Put each ufunc's samples in a column-wise matrix
        data = np.vstack([param._mcpts for param in params]).T

        # Apply the correlation matrix to the sampled data
        new_data = self.induce_correlations_SVD(data, correlation_matrix)

        # Re-set the samples to the respective variables
        for i in range(len(params)):
            params[i]._mcpts = new_data[:, i]

    def impl_SVD(self,A):
        """
        Calculate the SVD equivalent of the Cholesky lower triangular matrix.
        """
        A = np.array(A)
        assert A.shape[0] == A.shape[1], "Input matrix must be square"

        U, D, _ = np.linalg.svd(A)
        D = np.sqrt(np.eye(len(D))*D)
        L = np.dot(U,D)
        return np.array(L)
    
    def induce_correlations_SVD(self,data,correlation_matrix):
        # This method was found at https://github.com/tisimst/mcerp/blob/master/mcerp/correlate.py
        # only slightly modified to accomodate SVD when Cholesky decomposition fails with negative eigenvalues.
        """
        Induce a set of correlations on a column-wise dataset
    
        Parameters
        ----------
        data : 2d-array
            An m-by-n array where m is the number of samples and n is the
            number of independent variables, each column of the array corresponding
            to each variable
        correlation_matrix : 2d-array
            An n-by-n array that defines the desired correlation coefficients
            (between -1 and 1). Note: the matrix must be symmetric and
            positive-definite in order to induce.
    
        Returns
        -------
        new_data : 2d-array
            An m-by-n array that has the desired correlations.
        """
        # Create an rank-matrix
        data_rank = np.vstack([mc.rankdata(datai) for datai in data.T]).T

        # Generate van der Waerden scores
        data_rank_score = data_rank / (data_rank.shape[0] + 1.0)
        data_rank_score = norm(0, 1).ppf(data_rank_score)

        # Calculate the lower triangular matrix of the Cholesky decomposition
        # of the desired correlation matrix
        p = self.impl_SVD(correlation_matrix)

        # Calculate the current correlations
        t = np.corrcoef(data_rank_score, rowvar=0)

        # Calculate the lower triangular matrix of the Cholesky decomposition
        # of the current correlation matrix
        q = self.impl_SVD(t)

        # Calculate the re-correlation matrix
        s = np.dot(p, np.linalg.inv(q))

        # Calculate the re-sampled matrix
        new_data = np.dot(data_rank_score, s.T)

        # Create the new rank matrix
        new_data_rank = np.vstack([mc.rankdata(datai) for datai in new_data.T]).T

        # Sort the original data according to new_data_rank
        for i in range(data.shape[1]):
            _, order = np.unique(
                np.hstack((data_rank[:, i], new_data_rank[:, i])), return_inverse=True
            )
            old_order = order[: new_data_rank.shape[0]]
            new_order = order[-new_data_rank.shape[0] :]
            tmp = data[np.argsort(old_order), i][new_order]
            data[:, i] = tmp[:]

        return data

    def normal_lhs_monte_carlo(self): 
        """Returns non-standardized correlated normal distributions as columns in a pandas data frame. 
        The distributions are not truncated If another distribution shape is needed, 
        truncate and/or transform the distributions after applying this method. 
        kwarg: corr_matrix should be upper triangular in form of a list.
        E.g. matrix = [[NaN, 1, 1, 0],[1, NaN, 1, 0],[1, 1, NaN, 0],[0, 0, 0, NaN]]
         --> corr_matrix = [1,1,0,1,0,0]"""  

        input_quantities_dictionary = self.input_dictionary['input_quantities_dictionary']
        stdevs = self.partial_analysis_results['input_standard_uncertainty']
        draws = self.input_dictionary['monte_carlo_draws']

        if 'list_correlation_matrix' in self.input_dictionary:
            corr_matrix = self.input_dictionary['list_correlation_matrix']
        else: corr_matrix = None

        labels = list(input_quantities_dictionary.keys())
        correlated_distributions = list()
        applied_corr_matrices = list()

        
        mc.npts = draws

        for t in range(len(list(input_quantities_dictionary.values())[0])): 
            # Set up correlation matrix.
            if corr_matrix == None:
                c_matrix = np.eye(len(labels))
            else:
                c_matrix = self.expand_corr_matrix(labels,corr_matrix[t])
                

            means = list() 
            for a in range(len(labels)):
                val = list(input_quantities_dictionary.values())[a][t][0]
                means.append(val)
            distributions = list()
            mc_sims = list()
            for counter in range(len(labels)): 
                # Set up normal simulation model for each input variable and collect draws in variable to distributions list.
                input_var = mc.N(0,1) # Alternative ways to get scipy distributions: https://pythonhosted.org/mcerp/distribution_constructors.html
                mc_sims.append(input_var)
            # By definition all correlation matrices must be at least positive semidefinite, but small errors could throw it.
            # If it isn't also positive definite, the mcerp Cholesky algorithm won't work. 
            if not np.all(np.linalg.eigvals(c_matrix) > 0):
                print('The correlation matrix is not positive definite. Switching from Cholesky decomposition to SVD.')
                self.correlate_SVD(mc_sims,c_matrix)
                applied_corr_matrix = mc.correlation_matrix(mc_sims)
            else:
                # Impose correlation using Iman-Conover algorithm.
                mc.correlate(mc_sims,c_matrix)
                applied_corr_matrix = mc.correlation_matrix(mc_sims)

            # Transform the output to have correct mean and standard deviation, using means and stdevs and put into data frame 
            for counter in range(len(labels)):
                dist = mc_sims[counter]._mcpts
                new_dist = means[counter] + (dist - dist.mean())*(stdevs[counter][t]/dist.std())
                distributions.append(new_dist)

            correlated_distributions_t = pd.DataFrame(distributions).transpose()
            correlated_distributions_t.columns = labels
            correlated_distributions.append(correlated_distributions_t)
            applied_corr_matrices.append(applied_corr_matrix)
            
        return correlated_distributions,applied_corr_matrices 
    
    def compute_exprs(self,expression,input_quantities): 
        """ input arguments:
                - expression <str>: python valid math expression to be evaluated.
                - input_quantities <dict>: dictionary where keys are the input variable labels, and the entries are the corresponding nominal input values."""
        ALLOWED_NAMES = {k: v for k, v in math.__dict__.items() if not k.startswith("__")}
        ALLOWED_NAMES.update(input_quantities)
        code = compile(expression, "<string>", "eval")
        for name in code.co_names:
            if name not in ALLOWED_NAMES:
                raise NameError(f"The use of '{name}' is not allowed")
        return eval(code, {"__builtins__": {}}, ALLOWED_NAMES)

    def calculate_standard_uncertainty(self): 
        """The code assumes a coverage factor k = 2 and symmetric probability distributions.
        Allowed distribution is currently 'normal'. This may be updated in the future."""

        input_quantities_dictionary = self.input_dictionary['input_quantities_dictionary']
        probability_distributions = self.input_dictionary['probability_distributions']
        expanded_unc_dictionary = self.input_dictionary['expanded_unc_dictionary']

        expanded_uncertainty = list()
        var_no = 0
        for key in expanded_unc_dictionary.keys():
            content = expanded_unc_dictionary[key]
            exp_unc = list()
            for t in range(len(content)):
                exp_unc.append(content[t][0])
            var_no = var_no + 1
            expanded_uncertainty.append(exp_unc)
        
        a = 0
        t = 0
        labels = list(input_quantities_dictionary.keys())
        standard_uncertainty = list()
        for i in probability_distributions:
            std_uncertainty = list()
            for t in range(len(input_quantities_dictionary[labels[a]])):
                if i.lower() == "normal":
                    std_uncertainty.append((expanded_uncertainty[a][t]/2)) 
                else:
                    print(f"A {i.lower()} distribution for the variable {labels[a]} is not supported, assuming a normal distribution.") 
                    std_uncertainty.append((expanded_uncertainty[a][t]/2)) 
            standard_uncertainty.append(std_uncertainty)
            a = a + 1

        self.partial_analysis_results['input_standard_uncertainty'] = standard_uncertainty
        self.partial_analysis_results['input_expanded_uncertainty'] = expanded_uncertainty
        
    def calculate_sensitivity_coeff_and_correlation_term(self,labels,partial_derivatives, standard_uncertainties,corr_matrix=None): 
        
        correlation_term = list()
        sensitivity_coeff = [a*b for a,b in zip(partial_derivatives,standard_uncertainties)]
        if corr_matrix != None:
            counter = 0 
            for i in range(len(labels)-1):
                for j in range(len(labels)-1):
                       if j >= i:
                        correlation_term.append(
                            2*corr_matrix[counter]*partial_derivatives[i]*partial_derivatives[j+1]*standard_uncertainties[i]*standard_uncertainties[j+1])
                        counter = counter + 1
        return sensitivity_coeff, correlation_term

    def calculate_standard_combined_uncertainty(self,coeff,correlation_term): 
        square = [a*b for a,b in zip(coeff,coeff)]
        std_comb_u = math.sqrt(sum(square) + sum(correlation_term))
        return square, std_comb_u

    def expand_corr_matrix(self,labels,corr_matrix):
        tab = [[0]*len(labels) for _ in range(len(labels))]
        count = 0
        for i in range(len(labels)):
            for j in range(len(labels)):
                if i == j: tab[i][j] = 1
                elif j > i: # Filling in upper triangle of matrix.
                    tab[i][j] = corr_matrix[count]
                    count = count + 1
        i_lower = np.tril_indices(len(labels), -1)
        corr_matrix = np.array(tab)
        corr_matrix[i_lower]=corr_matrix.T[i_lower]
        return corr_matrix

    def perform_uncertainty_analysis(self):
        computation_method = self.input_dictionary['computation_method']

        if computation_method == "analytical":
            self.analytical()
        elif computation_method == "numerical":
            self.numerical()
        elif computation_method == 'lhs_montecarlo': 
            montecarlo_vectors, out_correlation_matrices = self.normal_lhs_monte_carlo() 

        if 'montecarlo' in computation_method:
            self.store_montecarlo_analysis_results(
                montecarlo_vectors, out_correlation_matrices)
        else:
            self.store_analytical_or_numerical_results()


    def store_analytical_or_numerical_results(self):

        inputs_df = self.add_input_sheet_time_series()
        outputs_df = self.add_output_sheet_time_series()

        # Getting necessary information from the uncertainty calculator object.
        input_quantities_dictionary = self.input_dictionary['input_quantities_dictionary']
        out_var_label = self.input_dictionary['out_var_label']
        functional_relationship = self.input_dictionary['functional_relationship']
        standard_uncertainties = self.partial_analysis_results['input_standard_uncertainty']
        if 'time_stamps' in self.input_dictionary.keys(): time_stamps = self.input_dictionary['time_stamps']
        else: time_stamps = pd.DataFrame()

        input_quantities_dictionary = self.input_dictionary['input_quantities_dictionary']
        out_var_label = self.input_dictionary['out_var_label']
        functional_relationship = self.input_dictionary['functional_relationship']
        standard_uncertainties = self.partial_analysis_results['input_standard_uncertainty']
        corr_matrix = self.input_dictionary['list_correlation_matrix']
        partial_derivatives = self.partial_analysis_results['partial_derivatives']
        partial_expressions = self.partial_analysis_results['partial_expressions']

        keys = list(input_quantities_dictionary.keys())
        input_units = list()
        for var in list(input_quantities_dictionary.values()):
            input_units.append(var[0][1])

        # Initializations. 
        standard_uncertainty_df = pd.DataFrame(columns=keys)
        formula_for_sensitivity_coeff_df = pd.DataFrame(columns=keys)
        sensitivity_coeff_df = pd.DataFrame(columns=keys)
        squared_contribution_df = pd.DataFrame(columns=keys)
        values_df = pd.DataFrame()
        standard_combined_uncertainties_df = pd.DataFrame()
        standard_rel_combined_uncertainties_df = pd.DataFrame()
        
        for col_idx in range(len(keys)):
            standard_uncertainty_df[keys[col_idx]] = standard_uncertainties[col_idx]
            formula_for_sensitivity_coeff_df[keys[col_idx]] = partial_expressions[col_idx]
            sens_coeffs = list()
            squares = list()
            values = list()
            std_combined_uncertainties = list()
            relative_std_combined_uncertainty = list()
            for t in range(len(list(input_quantities_dictionary.values())[0])):
                part_ders = list()
                std_unc = list()
                input_values = list()
                for a in range(len(keys)):
                    part_ders.append(partial_derivatives[a][t])
                    std_unc.append(standard_uncertainties[a][t])
                    input_values.append(list(input_quantities_dictionary.values())[a][t][0])
                coeff, correlation_term = self.calculate_sensitivity_coeff_and_correlation_term(
                     keys,part_ders, std_unc,corr_matrix=corr_matrix[t])
                input_quantities = dict(zip(keys,input_values))
                square, std_combined_uncertainty = self.calculate_standard_combined_uncertainty(
                    coeff,correlation_term)
                squares.append(square[col_idx])
                sens_coeffs.append(part_ders[col_idx])
                std_combined_uncertainties.append(std_combined_uncertainty)
                output_value = self.compute_exprs(functional_relationship,input_quantities)
                values.append(output_value)
                relative_std_combined_uncertainty.append(std_combined_uncertainty/output_value)
                
            squared_contribution_df[keys[col_idx]] = squares
            sensitivity_coeff_df[keys[col_idx]] = sens_coeffs
        standard_combined_uncertainties_df[out_var_label] = std_combined_uncertainties
        standard_rel_combined_uncertainties_df[out_var_label] = relative_std_combined_uncertainty
        values_df[out_var_label] = values

        units_df = pd.DataFrame(self.input_dictionary['out_var_unit'],columns=['units'], index=range(len(values_df)))
        absolute = pd.DataFrame({'Abs_'+out_var_label:len(values_df)*[True]})
        standard_uncertainty_df = pd.concat([time_stamps, standard_uncertainty_df], axis=1)
        formula_for_sensitivity_coeff_df = pd.concat(
            [time_stamps, formula_for_sensitivity_coeff_df], axis=1)
        sensitivity_coeff_df = pd.concat([time_stamps, sensitivity_coeff_df], axis=1)
        squared_contribution_df = pd.concat([time_stamps, squared_contribution_df], axis=1)
        standard_combined_uncertainties_df = pd.concat(
            [time_stamps, absolute, standard_combined_uncertainties_df,units_df], axis=1)
        standard_rel_combined_uncertainties_df = pd.concat(
            [time_stamps, standard_rel_combined_uncertainties_df], axis=1)
        values_df = pd.concat([time_stamps, values_df, units_df], axis=1)

        # Build out_correlation_matrix
        self.partial_analysis_results['correlation_matrix'] = corr_matrix
        correlation_matrix_df = self.add_correlations_sheet()
        
        # Store in uncertainty calculator object.
        self.analysis_results['output_values'] = values_df
        self.analysis_results['output_std_comb_uncertainties'] = standard_combined_uncertainties_df
        self.analysis_results['output_std_comb_rel_uncertainties'] = standard_rel_combined_uncertainties_df
        self.analysis_results['input_standard_uncertainties'] = standard_uncertainty_df
        self.analysis_results['formulae_for_sensitivity_coeff'] = formula_for_sensitivity_coeff_df
        self.analysis_results['sensitivity_coefficients'] = sensitivity_coeff_df
        self.analysis_results['squared_contribution'] = squared_contribution_df
        self.analysis_results['correlation_matrix'] = correlation_matrix_df

        if 'time_stamps' in self.input_dictionary:
            self.sheet_names = [
                'inputs','outputs','values','uncertainties','rel_uncertainties',
                'standard_input_uncertainty','formula_sens_coeff','sens_coeff','squared_contribution',
                'correlation_matrix']
            self.out_frames = [
                inputs_df, outputs_df, values_df,standard_combined_uncertainties_df,standard_rel_combined_uncertainties_df,
                standard_uncertainty_df,formula_for_sensitivity_coeff_df,sensitivity_coeff_df,squared_contribution_df,
                correlation_matrix_df]
        else: 
            sheet_names = ['inputs','outputs','correlation_matrix']
            outputs = self.add_output_sheet(
                values_df.iloc[0,0],values_df.columns[0],outputs_df.iloc[0,2],
                standard_combined_uncertainties_df.iloc[0,1],standard_rel_combined_uncertainties_df.iloc[0,0])
            computation_method = self.input_dictionary['computation_method']
            expanded_uncertainty = self.partial_analysis_results['input_expanded_uncertainty']
            
            partial_derivatives = [item for sublist in partial_derivatives for item in sublist]
            partial_expressions = [item for sublist in partial_expressions for item in sublist]
            expanded_uncertainty = [item for sublist in expanded_uncertainty for item in sublist]
            inputs = self.add_input_sheet(
                input_quantities_dictionary,expanded_uncertainty,inputs_df['ProbDist'].to_list(),standard_uncertainty_df.iloc[0,:].to_list(),
                computation_method, partial_derivatives=partial_derivatives, partial_expressions=partial_expressions,
                correlation_term = correlation_term, square=square)
            self.out_frames = [inputs, outputs, correlation_matrix_df]
            self.sheet_names = sheet_names

    def montecarlo_sensitivity_coeffs(self,montecarlos,draws=None):
        if not draws is None:
            self.montecarlo_sensitivity_draws = draws

        functional_relationship = self.input_dictionary['functional_relationship']

        labels = list()
        for elm in self.input_dictionary['labels_and_units']:
            labels.append(elm.split('|')[0])
        sensitivity_coeff_df = pd.DataFrame(columns=labels,index=range(len(montecarlos)))

        # montecarlos is a list of data frames, one list per time point.
        nominal_values = pd.DataFrame(columns=labels) 
        standard_deviations = pd.DataFrame(columns=labels, index = ['OAT','inp'])
        stdev_l = list()
        for t in range(len(montecarlos)):
            montecarlos_t = montecarlos[t][montecarlos[t].columns[:-1]].iloc[:self.montecarlo_sensitivity_draws,:]
            nominal_list = list()
            for var in range(len(montecarlos_t.columns)):
                nominal_list.append(self.input_dictionary['input_quantities_dictionary'][labels[var]][t][0])
            nominal_values.loc[0,:] = nominal_list
            nominal_values = nominal_values.loc[np.repeat(nominal_values.index, self.montecarlo_sensitivity_draws)].reset_index(drop=True)
            # We now have a data frame that contains nominal values montecarlo_sensitivity_draws times.
            # For each variable, switch it out one at a time with the corresponding variable from montecarlos.
            for var in montecarlos_t.columns.to_list():
                calculation_df = nominal_values.copy()
                calculation_df[var] = montecarlos_t[var]
                # Now compute the OAT stdev, and var stdev and find the sensitivity coefficient.
                # For each row in montecarlo_vectors calculate output.
                for row_n in range(len(calculation_df)):
                    vals = calculation_df.loc[row_n].tolist()
                    input_quantities = dict(zip(montecarlos_t.columns.to_list(),vals))
                    calculation_df.at[row_n,self.input_dictionary['out_var_label']] = self.compute_exprs(
                        functional_relationship,input_quantities)
        
                standard_deviations.loc['OAT',var] = calculation_df[self.input_dictionary['out_var_label']].std()
                standard_deviations.loc['inp',var] = calculation_df[var].std()
            sensitivity_coeff = (standard_deviations.loc['OAT']/standard_deviations.loc['inp']).tolist()
            stdev_l.append(standard_deviations)
            sensitivity_coeff_df.loc[t,:] = sensitivity_coeff
        
        if 'time_stamps' in self.input_dictionary:
            sensitivity_coeff_df = pd.concat((self.input_dictionary['time_stamps'],sensitivity_coeff_df),axis=1)
        return sensitivity_coeff_df

    # This function isn't entirely finished, only works for input variables with non-truncated normal distributions.
    def store_montecarlo_analysis_results(self,montecarlo_vectors,out_correlation_matrix):  
        """ So far this method assumes that the probability distributions of the input variables are normal. If this is not the case,
        the montecarlo draws in montecarlo_vectors should be transformed first."""
        
        inputs_df = self.add_input_sheet_time_series()
        outputs_df = self.add_output_sheet_time_series()
        
        # Initializing lists that will go to output data frames.
        values = list()
        std_combined_uncertainty = list()
        standard_rel_uncertainty = list()
        coverage_interval = list()
        standard_inp_unc = list()

        # Getting necessary information from the uncertainty calculator object.
        input_quantities_dictionary = self.input_dictionary['input_quantities_dictionary']
        out_var_label = self.input_dictionary['out_var_label']
        functional_relationship = self.input_dictionary['functional_relationship']
        standard_uncertainties = self.partial_analysis_results['input_standard_uncertainty']
        if 'time_stamps' in self.input_dictionary.keys(): time_stamps = self.input_dictionary['time_stamps']
        else: time_stamps = pd.DataFrame()

        t = 0
        for frame_t in montecarlo_vectors:
            frame_t[out_var_label] = ""
            std_uncertainties = list()
            # For each row in montecarlo_vectors calculate output.
            for row_n in range(len(frame_t)):
                vals = frame_t.loc[row_n].tolist()[:-1]
                input_quantities = dict(zip(list(input_quantities_dictionary.keys()),vals))
                frame_t.at[row_n,out_var_label] = self.compute_exprs(
                    functional_relationship,input_quantities)
            for a in range(len(input_quantities_dictionary.keys())):
                std_uncertainties.append(standard_uncertainties[a][t])

            # This needs to go to the output data frame.
            output_value = frame_t[out_var_label].mean()
            std_unc = frame_t[out_var_label].std()
            values.append(output_value)
            std_combined_uncertainty.append(std_unc)
            standard_rel_uncertainty.append(std_unc/output_value) 
            quantiles = frame_t[out_var_label].quantile([0.025,0.975])
            coverage_interval.append(
                str(output_value-quantiles.iloc[0])+' - '+str(output_value+quantiles.iloc[1]))
            standard_inp_unc.append(std_uncertainties)

            t = t + 1

        values_df = pd.DataFrame({out_var_label:values})
        standard_combined_uncertainties_df = pd.DataFrame({out_var_label:std_combined_uncertainty})
        standard_rel_combined_uncertainties_df = pd.DataFrame({out_var_label:standard_rel_uncertainty})
        coverage_interval_df = pd.DataFrame({out_var_label:coverage_interval})
        standard_uncertainty_df = pd.DataFrame(standard_inp_unc, columns=input_quantities_dictionary.keys())

        # Add time stamps and inform that the outputted uncertainty is absolute.
        units_df = pd.DataFrame(self.input_dictionary['out_var_unit'],columns=['units'], index=range(len(values_df)))
        absolute = pd.DataFrame({'Abs_'+out_var_label:len(values_df)*[True]})
        standard_uncertainty_df = pd.concat([time_stamps, standard_uncertainty_df], axis=1)
        coverage_interval_df = pd.concat([time_stamps,coverage_interval_df],axis=1)
        standard_combined_uncertainties_df = pd.concat([time_stamps, absolute, standard_combined_uncertainties_df,units_df], axis=1)
        standard_rel_combined_uncertainties_df = pd.concat([time_stamps, standard_rel_combined_uncertainties_df], axis=1)
        values_df = pd.concat([time_stamps, values_df, units_df], axis=1)

        # Correlation matrix
        matrix_list = list()
        for t in range(len(out_correlation_matrix)):
            vals = list()
            for i in range(len(out_correlation_matrix[t]-1)):
                vals.append(out_correlation_matrix[t][i][i+1:].tolist())
            matrix_list.append([item for sublist in vals for item in sublist]) 
        self.partial_analysis_results['correlation_matrix'] = matrix_list
        correlation_matrix_df = self.add_correlations_sheet()
        sensitivity_coefficients_df = self.montecarlo_sensitivity_coeffs(montecarlo_vectors)

        # Store analysis results to the uncertainty calculator object.
        # Store in uncertainty calculator object.
        self.analysis_results['output_values'] = values_df
        self.analysis_results['output_std_comb_uncertainties'] = standard_combined_uncertainties_df
        self.analysis_results['output_std_comb_rel_uncertainties'] = standard_rel_combined_uncertainties_df
        self.analysis_results['input_standard_uncertainties'] = standard_uncertainty_df
        self.analysis_results['coverage_interval'] = coverage_interval_df
        self.analysis_results['correlation_matrix'] = correlation_matrix_df
        self.analysis_results['sensitivity_coefficients'] = sensitivity_coefficients_df

        # Return data frames for excel output.
        if 'time_stamps' in self.input_dictionary.keys():
            self.sheet_names = [
                'inputs','outputs','values','uncertainties','rel_uncertainties','sensitivity_coefficients','coverage_interval','standard_uncertainty_input',
                'correlation_matrix']
            self.out_frames = [
                inputs_df, outputs_df, values_df, standard_combined_uncertainties_df,standard_rel_combined_uncertainties_df,
                sensitivity_coefficients_df, coverage_interval_df,standard_uncertainty_df,correlation_matrix_df]
        else: 
            self.sheet_names = ['inputs','outputs','correlation_matrix']
            computation_method = self.input_dictionary['computation_method']
            outputs = self.add_output_sheet(
                values_df.iloc[0,0],values_df.columns[0],outputs_df.iloc[0,2],
                standard_combined_uncertainties_df.iloc[0,1],standard_rel_combined_uncertainties_df.iloc[0,0],coverage_interval_df.iloc[0,0])
            expanded_uncertainty = self.partial_analysis_results['input_expanded_uncertainty']
            expanded_uncertainty = [item for sublist in expanded_uncertainty for item in sublist]
            inputs = self.add_input_sheet(
                input_quantities_dictionary,expanded_uncertainty, inputs_df['ProbDist'].to_list(),standard_uncertainty_df.iloc[0,:].to_list(),
                computation_method,sensitivity_coefficients_mc=sensitivity_coefficients_df)
            self.out_frames = [inputs, outputs, correlation_matrix_df ]

    def analytical(self): 

        input_quantities_dictionary = self.input_dictionary['input_quantities_dictionary']
        functional_relationship = self.input_dictionary['functional_relationship']

        partial_derivatives = list()
        partial_expressions = list()
        
        keys = list(input_quantities_dictionary.keys())
        # Make sure the input dictionary is in the form that compute_expr takes.
        
        for a in range(len(keys)):
            partial_ders = list()
            partial_exprs = list()
            for t in range(len(list(input_quantities_dictionary.values())[0])):
                vals = list()
                for b in range(len(keys)):   
                    vals.append(list(input_quantities_dictionary.values())[b][t][0])
                
                input_quantities = dict(zip(keys,vals))
            
                partialExp = str(diff(functional_relationship,keys[a]))
                partial_ders.append(self.compute_exprs(partialExp,input_quantities))
                partial_exprs.append(partialExp)
            
            partial_derivatives.append(partial_ders)
            partial_expressions.append(partial_exprs)
        
        self.partial_analysis_results['partial_derivatives'] = partial_derivatives
        self.partial_analysis_results['partial_expressions'] = partial_expressions
        
    def numerical(self,input_quantities_dictionary,functional_relationship,standard_uncertainties):

        input_quantities_dictionary = self.input_dictionary['input_quantities_dictionary']
        functional_relationship = self.input_dictionary['functional_relationship']
        standard_uncertainties = self.partial_analysis_results['input_standard_uncertainties']

        partial_derivatives = list()
        partial_expressions = list()
        
        keys = list(input_quantities_dictionary.keys())
        for a in range(len(keys)):
            partial_ders = list()
            partial_exprs = list()
            for t in range(len(list(input_quantities_dictionary.values())[0])):
                vals = list()
                for b in range(len(keys)):   
                    vals.append(list(input_quantities_dictionary.values())[b][t][0])
                
                input_quantities = dict(zip(keys,vals))

                add_unc = functional_relationship.replace(keys[a],"("+keys[a]+"+"+str(standard_uncertainties[a][t])+")")
                sub_unc = functional_relationship.replace(keys[a],"("+keys[a]+"-"+str(standard_uncertainties[a][t])+")")
                partialExp = "(("+add_unc+")-("+sub_unc+"))/(2*"+str(standard_uncertainties[a][t])+")"
                partial_ders.append(self.compute_exprs(partialExp,input_quantities))
                partial_exprs.append(partialExp)
            partial_derivatives.append(partial_ders)
            partial_expressions.append(partial_exprs)
        self.partial_analysis_results['partial_derivatives'] = partial_derivatives
        self.partial_analysis_results['partial_expressions'] = partial_expressions

    def output_to_excel(self, path, filename): 
        data_frames = self.out_frames  
        sheet_names = self.sheet_names
        
        computation_method = self.input_dictionary['computation_method']
        functional_relationship = self.input_dictionary['functional_relationship']
        monte_carlo_draws = self.input_dictionary['monte_carlo_draws']

        # Inform where the output can be found.
        replace_strings = ['_input_','_input','input_','input',computation_method,'.xlsx']
        for replace_string in replace_strings:
            if replace_string in filename:
                filename = filename.replace(replace_string,'')
        
        filename = path+ '/output_'+filename+'_'+str(computation_method)
        print("Writing results to "+filename+".xlsx")
        
        # Add settings
        if any(self.partial_analysis_results['correlations_exist']) == True:
            corrs = 'yes'
        else: corrs = 'no'

        
        settings_df = pd.DataFrame() 
        settings_df["Inputs"] = [
            "Functional relationship:", "Method of computation:", 
            "There are correlations:"]
        settings_df["Options"] = [functional_relationship, computation_method, corrs]
        if 'montecarlo' in computation_method:
            settings_df.loc[len(settings_df.index)] = ["Monte Carlo draws: ", monte_carlo_draws]

        # write data frames to excel
        correlation_matrix_df = data_frames[-1]
        data_frames = data_frames[:-1]
        with pd.ExcelWriter(filename+".xlsx") as writer:
            settings_df.to_excel(writer, sheet_name='settings',index=False)  
            i = 0
            for frame in data_frames:
                frame.to_excel(writer, sheet_name=sheet_names[i],index=False)
                i = i + 1
            if corrs == 'yes':
                correlation_matrix_df.to_excel(writer, sheet_name='correlation_matrix',index=True)
            else: correlation_matrix_df.to_excel(writer, sheet_name='correlation_matrix',index=False)

    def add_input_sheet_time_series(self):
        input_quantities_dictionary = self.input_dictionary['input_quantities_dictionary']
        keys = list(input_quantities_dictionary.keys())
        input_units = list()
        for var in list(input_quantities_dictionary.values()):
            input_units.append(var[0][1])

        inputs_df = pd.DataFrame() 
        inputs_df["Label"] = keys
        inputs_df["Unit"] = input_units
        inputs_df["ProbDist"] = self.input_dictionary['probability_distributions']
        return inputs_df

    def add_output_sheet_time_series(self):
        label = self.input_dictionary['out_var_label']
        unit = self.input_dictionary['out_var_unit']
        outputs_df = pd.DataFrame() 
        outputs_df["Outputs"] = ["Output variable:","Std. comb. unc:","Std. rel. comb. unc."]
        outputs_df["Label"] = [label,'U_{'+label+'}','u_rel_{'+label+'}']
        outputs_df["Unit"] = [unit,unit,'']
        return outputs_df

    def add_input_sheet(
        self,input_quantities_dictionary,expanded_uncertainty,probability_distributions,standard_uncertainties,
        computation_method,partial_derivatives=None,partial_expressions=None,
        correlation_term=None,square=None,sensitivity_coefficients_mc=None):
        labels = list(input_quantities_dictionary.keys())
        values = list()
        units = list()
        for elm in list(input_quantities_dictionary.values()):
            values.append(elm[0][0])
            units.append(elm[0][1])
        
        inputs_df = pd.DataFrame()
        
        inputs_df["Input parameters \n X"] = [m+" ["+n+"]" for m,n in zip(labels,units)]
        inputs_df["Input values \n x_i"] = values
        inputs_df["Expanded rel. uncertainty \n U(x_i)/x_i"] = expanded_uncertainty
        inputs_df["Distribtution"] = probability_distributions
        inputs_df["Standard uncertainty \n u(x_i)"] = standard_uncertainties
        if not "montecarlo" in computation_method:
            inputs_df["Sensitivity coeff. \n c_i"] = partial_derivatives
            inputs_df["Formula for c_i = (df/dx_i)"] = partial_expressions
            if self.partial_analysis_results['correlations_exist'][0] == True:
                inputs_df["(c_i*u(v_i))^2 + 2*R(x_i,x_j)*(df/dx_i)*(df/dx_j)*u(x_i)*u(x_j)"] = [
                    m+n for m,n in zip(square,correlation_term)]
            else:
                inputs_df["(c_i*u(v_i))^2 + 2*R(x_i,x_j)*(df/dx_i)*(df/dx_j)*u(x_i)*u(x_j)"] = square  
        else: inputs_df["Sensitivity coeff."] = sensitivity_coefficients_mc.iloc[0,:].to_list() 
        return inputs_df
        
    def add_correlations_sheet(self):
        input_labels = list(self.input_dictionary['input_quantities_dictionary'].keys())
        out_correlation_matrix = self.partial_analysis_results['correlation_matrix']
        if  'time_stamps' in self.input_dictionary.keys():
            time_stamps = self.input_dictionary['time_stamps']

        if not 'time_stamps' in self.input_dictionary.keys(): 
            if self.partial_analysis_results['correlations_exist'][0] == False:
                correlation_matrix_df = pd.DataFrame()
                correlation_matrix_df["Correlation matrix"] = ["No correlation in input variables."]
            else:
                correlation_matrix_df = pd.DataFrame(self.expand_corr_matrix(input_labels, out_correlation_matrix[0]))
                correlation_matrix_df.columns = input_labels
                correlation_matrix_df.index = input_labels
        else:
            correlation_matrix_df = pd.DataFrame(out_correlation_matrix)
            col_list = list()
            combi = list()
            for label in input_labels:
                for next_label in input_labels:
                    if not label == next_label and [next_label,label] not in combi:
                        col_list.append('R('+label+','+next_label+')')
                        combi.append([label,next_label])
            correlation_matrix_df.columns = col_list
            correlation_matrix_df = pd.concat([time_stamps, correlation_matrix_df],axis=1)
        return correlation_matrix_df
        
    def add_output_sheet(
        self,output_value,output_label,output_unit,std_combined_uncertainty,
        relative_std_combined_uncertainty,coverage_interval=None):
        outputs = pd.DataFrame() 
        outputs["Outputs"] = ["Output value:", "Std. comb. unc.:", "Rel. std. comb. unc.:"]
        outputs["Label"] = [output_label,'U_{'+output_label+'}','u_rel_{'+output_label+'}']
        outputs["Value"] = [output_value,std_combined_uncertainty,relative_std_combined_uncertainty]
        outputs["Unit"] = [output_unit, output_unit,"-"]
        if coverage_interval != None:
            outputs.loc[len(outputs.index)] = ["{:.0%} coverage interval".format(0.95), '', coverage_interval, output_unit]
        return outputs     

    def get_inputs_from_file(self,input_path,input_file,time_series=False): 
        """It is not necessary to read inputs from a specific file format to run the calculator, but the input information must be passed
        to the calculator by means of a dictionary that is similar to the one returner from this method."""

        
        xls = pd.ExcelFile(input_path  + "/"+ input_file)

        # Get inputs
        warnings.simplefilter("ignore")
        df_inputs = pd.read_excel(xls, 'inputs')
        warnings.resetwarnings()
        input_quantity_labels = df_inputs["Label"].tolist()
        units = df_inputs["Unit"].tolist()
        probability_distributions = df_inputs["ProbDist"].tolist()
        
        
        if time_series == False:
            absoluteUncs = df_inputs["AbsUnc"].to_list()
            values = df_inputs["Value"].tolist()
            expanded_uncertainty = df_inputs["ExpandedUnc"].tolist()
            for idx in range(len(expanded_uncertainty)):
                if absoluteUncs[idx] == False:
                    expanded_uncertainty[idx] = values[idx]*expanded_uncertainty[idx]
            expanded_unc_dictionary= dict(zip(input_quantity_labels,zip(expanded_uncertainty,units)))
            input_quantities_dictionary = dict(zip(input_quantity_labels,zip(values,units)))
            for key in input_quantities_dictionary.keys():
                expanded_unc_dictionary[key] = [list(expanded_unc_dictionary[key])]
                input_quantities_dictionary[key] = [list(input_quantities_dictionary[key])]
        else:
            values = list()
            expanded_uncertainty = list()
            df_values = pd.read_excel(xls,'values')
            df_unc = pd.read_excel(xls,'uncertainties')
            for label in input_quantity_labels:
                if 'Time_factor' in df_values: 
                    vals = [a*b for a,b in zip(df_values['Time_factor'].to_list(),df_values[label].to_list())]
                else: vals = df_values[label].to_list()
                values.append(vals)
                exp_unc = list()
                for t in range(len(df_unc[label])):
                    if df_unc['Abs_'+label].to_list()[t] == True:
                        exp_unc.append(df_unc['Time_factor'][t]*df_unc[label][t])
                    else: exp_unc.append(vals[t]*df_unc[label][t])
                expanded_uncertainty.append(exp_unc)
        input_quantities_dictionary = dict()
        expanded_unc_dictionary = dict()
        for t in range(len(values)):
            if isinstance(values[0],list):
                value_tuples = list(zip(values[t],len(values[0])*[units[t]]))
                exp_unc_tuples = list(zip(expanded_uncertainty[t],len(expanded_uncertainty[0])*[units[t]]))
                value_list = list()
                exp_unc_list = list()
                for elm in value_tuples:
                    value_list.append(list(elm))
                for elm in exp_unc_tuples:
                    exp_unc_list.append(list(elm))
                input_quantities_dictionary[input_quantity_labels[t]] = value_list
                expanded_unc_dictionary[input_quantity_labels[t]] = exp_unc_list
            else: 
                value_tuples = list(zip(values,units))
                exp_unc_tuples = list(zip(expanded_uncertainty,units))
                input_quantities_dictionary[input_quantity_labels[t]] = [list(value_tuples[t])]
                expanded_unc_dictionary[input_quantity_labels[t]] = [list(exp_unc_tuples[t])]
        

        # outputs
        df_outputs = pd.read_excel(xls, 'outputs')
        out_var_label = df_outputs["Label"].tolist()[0]
        out_var_unit = df_outputs["Unit"].tolist()[0]
        
        #settings
        warnings.simplefilter("ignore")
        df_settings = pd.read_excel(xls, 'settings')
        warnings.resetwarnings()
        functional_relationship = df_settings.iat[0,1]
        computation_method = df_settings.iat[1,1]
        monte_carlo_draws = df_settings.iat[2,1]
        

        # correlation matrix
        if time_series == False:
            list_correlation_matrix = list()
            df_correlation_matrix = pd.read_excel(xls, 'correlation_matrix')
            for column in df_correlation_matrix:
                idx = np.where(df_correlation_matrix[column].isna())[0]
                if len(idx) > 0:
                    lst = df_correlation_matrix[column].tolist()
                    list_correlation_matrix=list_correlation_matrix+lst[idx[0]+1:]
            list_correlation_matrix = [list_correlation_matrix]
        else:
            list_correlation_matrix = list()
            df_correlation_matrix = pd.read_excel(xls, 'correlation_matrix')
            for t in range(len(df_correlation_matrix)):
                list_correlation_matrix.append(df_correlation_matrix.iloc[t,1:].to_list())

        if time_series == True:
            time_stamps = pd.read_excel(xls,'values')
            time_stamps = time_stamps.iloc[:,:1]
            in_d = {
            'input_quantities_dictionary':input_quantities_dictionary,
            'expanded_unc_dictionary':expanded_unc_dictionary,'probability_distributions':probability_distributions,
            'out_var_label':out_var_label,'out_var_unit':out_var_unit,'functional_relationship':functional_relationship,
            'computation_method':computation_method,'monte_carlo_draws':monte_carlo_draws,
            'list_correlation_matrix':list_correlation_matrix,'time_stamps':time_stamps}
        else:
            in_d = {
            'input_quantities_dictionary':input_quantities_dictionary,
            'expanded_unc_dictionary':expanded_unc_dictionary,'probability_distributions':probability_distributions,
            'out_var_label':out_var_label,'out_var_unit':out_var_unit,'functional_relationship':functional_relationship,
            'computation_method':computation_method,'monte_carlo_draws':monte_carlo_draws,
            'list_correlation_matrix':list_correlation_matrix}
        
        self.input_dictionary = in_d
    
    def set_input_quantities(self,labels_and_units,values,time_stamps=None):
        """examples: labels_and_units = ['V1|m^3','V2|m^3','rho|(kg)/(m^3)']
                     values = [[300e3,301e3,301e3],[180e3,179e3,180e3],[810.0,810.0,810.0]]
                     
                     - Each list in values connects to the corresponding label and unit in labels_and_units.
                     - If time_stamps are not set, or does not match the number of samples in the input
                       variables, and the input data has more than one sample of each variable,
                       time_stamps is set to a default list ['t1','t2',...]. """
        
        input_quantities_dictionary = {}
        for elmno in range(len(labels_and_units)):
            elm = labels_and_units[elmno].split('|')
            input_quantities_dictionary[elm[0]] = []
            for t in range(len(values[0])):
                input_quantities_dictionary[elm[0]].append([values[elmno][t],elm[1]])

        if len(values[0])>1 and time_stamps is None:
            self.input_dictionary['time_stamps']=pd.DataFrame(['t'+str(elm) for elm in list(range(1,len(values[0])+1))],columns=['time_stamp'])
        if time_stamps is not None:
            if len(values[0])>1 and len(time_stamps)!=len(values[0]):
                self.input_dictionary['time_stamps']=pd.DataFrame(['t'+str(elm) for elm in list(range(1,len(values[0])+1))],columns=['time_stamp'])
            elif len(values[0])>1 and len(time_stamps)==len(values[0]):
                self.input_dictionary['time_stamps'] = time_stamps

        self.input_dictionary['input_quantities_dictionary'] = input_quantities_dictionary
        self.input_dictionary['labels_and_units'] = labels_and_units
    
    def set_input_expanded_uncertainties(self,values,absolute=True,distributions=None):
        """A coverage factor of k=2 is assumed.
        example input:  values = [[9000,9000,9000],[5400,5400,5400],[0.81,0.81,0.81]]

                     - Each list in values connects to the corresponding label and unit in self.labels_and_units.
                     - Indicate whether or not the inputted expanded uncertainties are absolute or not.
                     - If distributions=None, normal distributions are assumed, otherwise, indicate one type of distribution
                     for each input variable, f.ex distribution=['normal','normal','normal']. The length of distribution
                     must in that case equal the number of input variables."""

        allowed_distributions = ['normal']

        labels_and_units = self.input_dictionary['labels_and_units']
        expanded_unc_dictionary = {}
        for elmno in range(len(labels_and_units)):
            elm = labels_and_units[elmno].split('|')
            expanded_unc_dictionary[elm[0]] = []
            for t in range(len(values[0])):
                if absolute == True:
                    expanded_unc_dictionary[elm[0]].append([values[elmno][t],elm[1]])
                if absolute == False:
                    expanded_unc_dictionary[elm[0]].append(
                        [values[elmno][t]*self.input_dictionary['input_quantities_dictionary'][elm[0]][t][0],elm[1]])

        self.input_dictionary['probability_distributions'] = list()
        if  distributions is None:
            self.input_dictionary['probability_distributions']=len(labels_and_units)*['normal']
        else:
            if len(distributions)!=len(labels_and_units):
                self.input_dictionary['probability_distributions']=len(labels_and_units)*['normal']
            else: 
                for distno in range(len(distributions)):
                    if distributions[distno] in allowed_distributions:
                        self.input_dictionary['probability_distributions'].append(distributions[distno])
                    else: self.input_dictionary['probability_distributions'].append('normal')
        
        self.input_dictionary['expanded_unc_dictionary'] = expanded_unc_dictionary
    
    def set_functional_relationship(self,functional_relationship):
        """The functional relationship must contain all labels in self.input_dictionary['labels_and_units'].
        It must be a Python valid math expression. To avoid potential confusion, use plenty if parantheses."""

        ok = True
        labels = [elm.split("|")[0] for elm in self.input_dictionary['labels_and_units']]
        if not all(str(l) in functional_relationship for l in labels):
            ok = False
        if ok:
            self.input_dictionary['functional_relationship'] = functional_relationship
        else: raise Exception("Not all input labels are contained in the functional relationship.")
    
    def set_computational_method(self,computation_method=None, monte_carlo_draws=10000):
        """If no computational method is set, or an invalid computation method is set,
        the default method is 'lhs_montecarlo'.
        Other options are 'analytical' and 'numerical'. """
        if computation_method is None: computation_method = 'lhs_montecarlo'
        self.input_dictionary['computation_method'] = computation_method
        self.input_dictionary['monte_carlo_draws'] = monte_carlo_draws
      
    def set_output_labels_and_units(self,output_labels_and_units):
        """example: labels_and_units = ['M|kg']"""
        if isinstance(output_labels_and_units,list):
            self.input_dictionary['out_var_label'] =  output_labels_and_units[0].split('|')[0]
            self.input_dictionary['out_var_unit'] =  output_labels_and_units[0].split('|')[1]
        if isinstance(output_labels_and_units,str):
            self.input_dictionary['out_var_label'] =  output_labels_and_units.split('|')[0]
            self.input_dictionary['out_var_unit'] =  output_labels_and_units.split('|')[1]

    def set_correlation_matrix(self,correlation_matrix=None):
        """example: correlation_matrix=[[0.9,0,0],[0,0,0],[0,0,0]].
                    Each list corresponds to an upper triangular matrix, such that if 
                    self.input_dictionary['labels_and_units'] = ['V1|kg','V2|kg','rho|(kg)/(m^3)'], 
                    there is correlation between V1 and V2 at only at t1.
                    If ignored, zero correlation is assumed.
                    Zero correlation is also assumed if the correlation matrix does not 
                    match number of variables or length of time series"""
        corr_mat = list()
        time_len = len(self.input_dictionary['input_quantities_dictionary'][self.input_dictionary['labels_and_units'][0].split('|')[0]])
        if correlation_matrix is None or len(correlation_matrix)!=time_len:
            for _ in range(len(self.input_dictionary['labels_and_units'])):
                corr_mat.append([0]*time_len)
        if len(corr_mat)==0:
            for i in range(time_len):
                if len(correlation_matrix[i])!=len(self.input_dictionary['labels_and_units']):
                    for _ in range(len(self.input_dictionary['labels_and_units'])):
                        corr_mat.append([0]*time_len)
        if len(corr_mat)==0:corr_mat = correlation_matrix
        self.input_dictionary['list_correlation_matrix'] = corr_mat


