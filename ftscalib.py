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

def get_header(lines_in):
    """Extracts the header lines at the start of aln and cln files.
    Returns them in a single string"""    
    
    header = ""
    
    for header_line in lines_in[:4]:
        header += header_line
        
    return header
    

def get_lines(lines_in):
    """Assigns each line in an aln or cln file to a Line object. 
    Returns a list of Line objects"""
    
    all_lines = []
    
    for line in lines_in[4:]:
        extracted_line = Line(line)  # create object of class Line
        all_lines.append(extracted_line)
        
    return all_lines
        
        
def remove_bad_lines(line_list, snr_discrim):
    """Removes Line objects from line_list if they are not lsqfitted and/or
    have peak heights less than snr_discrim"""    
    
    good_lines = []
    
    for line in line_list:  # only keep if lsqfitted and snr > discrinator
        if line.tags == "L" and float(line.peak) >= snr_discrim:
            good_lines.append(line)
    return good_lines
    
    
def get_common_lines(aln_lines, cln_lines, wavenum_discrim):
    """Returns list of aln lines that match lines in the cln line_list within
    a tolerance of +- wavenum_discrim"""
    
    common_lines = []   
    
    for aln_line in aln_lines:
        no_of_matches = 0
        
        for cln_line in cln_lines:
            diff = float(aln_line.wavenumber) - float(cln_line.wavenumber)
            
            if abs(diff) <= wavenum_discrim:
                no_of_matches += 1
                delta_sig = diff / float(aln_line.wavenumber)  # delta sigma over sigma 
                common_lines.append((aln_line, delta_sig))
                
        if no_of_matches > 1:  # Raises exception if multiple possible matches found in cln file
            raise ValueError("Multiple line assignments for line: %s" % aln_line.wavenumber)
    
    return common_lines
    
                
def get_corr_factor(common_lines, std_dev_discrim):
    """Discards lines outside of the std_dev_discrim limits and recalculates
    correction factor, std dev and std error for the remaining lines. Returns
    lines used in the new calcuation, the rejected lines, newly calculated
    correction factor, std dev and std error."""
     
    lines_used = common_lines
    lines_rejected = []    
    corr_factor = calc_corr_factor(lines_used)
    std_dev, std_error = calc_std_dev(lines_used)
    limit = std_dev_discrim * std_dev
         
    for line in common_lines:           
        if line[1] > (corr_factor + limit) or line[1] < (corr_factor - limit):
            lines_rejected.append(line)
            lines_used.remove(line)
    
    corr_factor = calc_corr_factor(lines_used)
    print lines_used
    std_dev, std_error = calc_std_dev(lines_used)
     
    return lines_used, lines_rejected, corr_factor, std_dev, std_error    
    

def calc_corr_factor(common_lines):
    """Returns the weighted average of delta sigma over sigma for all good,
    common lines. This correction factor is weighted against SNR of the 
    lines"""
    
    weighted_factors = 0.   
    total_weight = 0.
    
    for line in common_lines:
        delta_sig = line[1]
        snr = float(line[0].peak)
        weighted_factors += (delta_sig * snr)
        total_weight += snr
        
    return weighted_factors / total_weight
    
def calc_corr_factor_sasha(common_lines):
    """Returns the weighted average of delta sigma over sigma for all good,
    common lines. This correction factor is weighted against SNR of the 
    lines"""
    
    weighted_factors = 0.   
    total_weight = 0.
    
    for line in common_lines:
        delta_sig = line[1]
     #   unc_aln = 
      #  unc_cln = 
        
        
        weighted_factors += (delta_sig * snr)
        total_weight += snr
        
    return weighted_factors / total_weight    
    
    
def calc_std_dev(common_lines):
    """Returns the standard deviation and standard error of the delta sigma 
    over sigma values of an output array from get_corr_factor"""
    
    delta_sigmas = []
    
    for line in common_lines:
        delta_sigmas.append(line[1])  # line[1] = delta sigma / sigma values 
        
    std_dev_array = np.array(delta_sigmas)  # convert to numpy array
    std_dev = np.std(std_dev_array)
    std_error = std_dev / np.sqrt(len(std_dev_array))
    
    return std_dev, std_error
    
    
#def calc_keff_unc_sasha():
    
    
    
    
#def calc_keff_unc_gill():
    
    
    
    

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

snr_discrim = 10.
wavenum_discrim = 0.05
std_dev_discrim = 3.5
lines = []

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

common_lines = get_common_lines(aln_good_lines, cln_good_lines, wavenum_discrim)

lines_used, lines_rejected, corr_factor, std_dev, std_error = get_corr_factor(common_lines, std_dev_discrim)
print len(lines_used), corr_factor, std_dev

aln_in.close()
cln_in.close()   
