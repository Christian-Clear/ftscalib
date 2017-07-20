# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 10:43:26 2017

@author: christian
"""

# -*- coding: utf-8 -*-
"""
Created on Thu May 25 14:20:06 2017

@author: christian


ftscalib

This program calibrates an xGremlin writelines output file against a previously
calibrated writelines file (or a set of standards - which must be in the 
writelines output format)
"""

import numpy as np

def get_header(header_lines_in):
    """Extracts the header lines at the start of aln and cln files.
    Returns them in a single string"""    
    
    header = ""
    
    for header_line in header_lines_in[:4]:
        header += header_line
        
    return header
    

def get_lines(get_lines_in):
    """Assigns each line in an aln or cln file to a Line object. 
    Returns a list of Line objects"""
    
    all_lines = []
    
    for line in get_lines_in[4:]:
        extracted_line = Line(line)  # create object of class Line
        all_lines.append(extracted_line)
        
    return all_lines
        
        
def remove_bad_lines(rm_bad_line_list, rm_bad_snr_disc):
    """Removes Line objects from line_list if they are not lsqfitted and/or
    have peak heights less than snr_discrim"""    
    
    good_lines = []
    
    for line in rm_bad_line_list:  # only keep if lsqfitted and snr > discrinator
        if line.tags == "L" and float(line.peak) >= rm_bad_snr_disc:
            good_lines.append(line)
    
    return good_lines
    
    
def get_common_lines(get_cmn_aln_lines, get_cmn_cln_lines, get_cmn_wavenum_discrim):
    """Returns list of aln lines that match lines in the cln line_list within
    a tolerance of +- wavenum_discrim"""
    
    get_cmn_common_lines = []   
    
    for aln_line in get_cmn_aln_lines:
        no_of_matches = 0
        
        for cln_line in get_cmn_cln_lines:
            diff = float(cln_line.wavenumber) - float(aln_line.wavenumber)
            
            if abs(diff) <= get_cmn_wavenum_discrim:
                no_of_matches += 1
                delta_sig = diff / float(aln_line.wavenumber)  # delta sigma over sigma 
                get_cmn_common_lines.append((aln_line, cln_line, delta_sig))
                
        if no_of_matches > 1:  # Raises exception if multiple possible matches found in cln file
            raise ValueError("Multiple line assignments for line: %s" % aln_line.wavenumber)
    
    return get_cmn_common_lines
    
                
 
def get_corr_factor_sasha(corr_fact_sasha_common_lines, std_dev_discrim):
    """Discards lines outside of the std_dev_discrim limits and recalculates
    correction factor, std dev and std error for the remaining lines. Returns
    lines used in the new calcuation, the rejected lines, newly calculated
    correction factor, std dev and std error."""
     
    lines_used = corr_fact_sasha_common_lines
    #print "lines used before:" + str(len(lines_used))
    lines_rejected = []    
    corr_factor, x = calc_corr_factor_sasha(lines_used)
    print "corr_factor_sasha: ", corr_factor
    std_dev, std_error = calc_std_dev(lines_used)
    limit = std_dev_discrim * std_dev
         
    for line in corr_fact_sasha_common_lines:   
        delta_sig = line[2]        
        if delta_sig > (corr_factor + limit) or delta_sig < (corr_factor - limit):
            lines_rejected.append(line)
            lines_used.remove(line)
    
    #print "lines used after:" + str(len(lines_used))    
    
    corr_factor, weights = calc_corr_factor_sasha(lines_used)
    corr_factor_unc = calc_corr_factor_unc_sasha(weights, corr_factor)
     
    return lines_used, lines_rejected, corr_factor, corr_factor_unc


def calc_corr_factor_sasha(calc_corr_factor_sasha_common_lines):
    """Returns the weighted average of delta sigma over sigma for all good,
    common lines. This correction factor is weighted against SNR of the 
    lines"""
    
    weights = []
    weighted_ki = 0.
    total_weight = 0.
    
    for line in calc_corr_factor_sasha_common_lines:
        ki = line[2]
        #print ki
        aln_line = line[0]
        cln_line = line[1]
        
        cln_unc = (float(cln_line.width) / 1000) / (2 * float(cln_line.peak))
        aln_unc = (float(aln_line.width) / 1000) / (2 * float(aln_line.peak))
        
        
        ki_unc = np.sqrt(np.power(aln_unc, 2) + np.power(cln_unc, 2))
        
        wi = np.power(ki_unc, -2)
        
        weighted_ki += (ki * wi)
        total_weight += wi
        

        weights.append((ki, wi))
           
    keff = (weighted_ki / total_weight)  # keff is the correction factor

    return keff, weights

#def calc_corr_factor_sasha(calc_corr_common_lines):
#    """Returns the weighted average of delta sigma over sigma for all good,
#    common lines. This correction factor is weighted against SNR of the 
#    lines"""
#    
#    weighted_factors = 0.   
#    total_weight = 0.
#    weights = []
#    
#    for line in calc_corr_common_lines:
#        ki = line[2]
#        aln_snr = float(line[0].peak)
#        aln_FWHM = float(line[0].width)
#        aln_unc = (aln_FWHM / (2 * aln_snr))
#        
#        cln_snr = float(line[1].peak)
#        cln_FWHM = float(line[1].width)
#        cln_unc = cln_FWHM / (2 * cln_snr) / float(line[1].wavenumber)
#        
#        tot_unc = np.sqrt(np.power(aln_unc, 2) + np.power(cln_unc, 2))
#        
#        wi = np.power(tot_unc, -2)
#        
#        weighted_factors += (ki * wi)
#        total_weight += wi
#        
#        weights.append((ki, wi))
#        
#    return (weighted_factors / total_weight), weights

    
    
def calc_corr_factor_unc_sasha(weights, keff):
    
    total_weight = 0. 
    keff_unc = 0.
    
    for weight in weights:
        ki = weight[0]
        wi = weight[1]        
            
        keff_unc += (wi + (np.power(wi, 2) * np.power(keff - ki, 2)))  # From Sahsa's paper 
        total_weight += wi
    
    keff_unc = (keff_unc**0.5) / total_weight
    
    return keff_unc   
    
    
def calc_std_dev(calc_std_dev_common_lines):
    """Returns the standard deviation and standard error of the delta sigma 
    over sigma values of an output array from get_corr_factor"""
    
    delta_sigmas = []
    
    for line in calc_std_dev_common_lines:
        delta_sigmas.append(line[2])  # line[1] = delta sigma / sigma values 
        
    std_dev_array = np.array(delta_sigmas)  # convert to numpy array
    std_dev = np.std(std_dev_array)
    std_error = std_dev / np.sqrt(len(std_dev_array))
    
    return std_dev, std_error
    
    
    

class Line:
    
    """Assigns all parameters from a line from an xgremlin "writelines" 
    command output file to the variables of a Line class object"""
    
    def __init__(self, xgremlin_line):  # assign values to variables from line string
    
        self.line = xgremlin_line[2:6]        
        self.wavenumber = xgremlin_line[8:20]        
        self.peak = xgremlin_line[21:30]        
        self.width = xgremlin_line[33:39]        
        self.dmp = xgremlin_line[42:48]        
        self.eq_width = xgremlin_line[49:59]        
        self.itn = xgremlin_line[64:65]        
        self.H = xgremlin_line[67:69]        
        self.tags = xgremlin_line[73:74]        
        self.epstot = xgremlin_line[75:85]        
        self.epsevn = xgremlin_line[86:96]        
        self.epsodd = xgremlin_line[97:107]        
        self.epsran = xgremlin_line[108:118]        
        self.id = xgremlin_line[119:145]        
        self.wavelength = xgremlin_line[150:160]
        
        

#### Main Program ####
def main():
        
    snr_discrim = 10.
    wavenum_discrim = 0.05
    std_dev_discrim = 0.5
    
    aln_in = open("test.aln", 'r')
    cln_in = open("test.cln", "r")
    
    aln_in_lines = aln_in.readlines()
    cln_in_lines = cln_in.readlines()
    
    aln_header = get_header(aln_in_lines)
    cln_header = get_header(cln_in_lines)
    
    aln_raw_lines = get_lines(aln_in_lines)
    cln_raw_lines = get_lines(cln_in_lines)
    
    aln_good_lines = remove_bad_lines(aln_raw_lines, snr_discrim)
    cln_good_lines = remove_bad_lines(cln_raw_lines, snr_discrim)
    
    cmn_lines = get_common_lines(aln_good_lines, cln_good_lines, wavenum_discrim)
    
   # common_lines2 = common_lines
    
    lines_used, lines_rejected, corr_factor, corr_factor_unc = get_corr_factor_sasha(cmn_lines, std_dev_discrim)
    print len(lines_used), corr_factor, corr_factor_unc
      
    aln_in.close()
    cln_in.close()   

main()