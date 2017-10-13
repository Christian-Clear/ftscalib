#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu May 25 14:20:06 2017

@author: Christian Clear


ftscalib

This program calibrates an xGremlin writelines output file against a previously
calibrated writelines file (or a set of standards - which must be in the 
writelines output format)
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os


    
def print_docstring():  # Main document string. Gets called if no arguements are passed to ftscalib
    docstring = """
    
    ========== ftscalib (v1.1) ============
    
    author : Christian Clear
    email  : cpc14@ic.ac.uk
    date   : 13/10/2017
    
    ftscalib is a python script for calibrating xGremlin .aln files against previously calibrated lines or Ar II standards.
    Use only the Ar II standards bundled with ftscalib or match the format exactly if using your own. ftscalib is designed 
    to be run in a unix terminal. 
    
    --------------------------------------------------------------------------------    
    Syntax : ftscalib <inp file>
    <inp file>  : The ftscalib (.inp) input file

    --------------------------------------------------------------------------------
    Outputs:
        calibrated linelist (.cln)
        uncertainties of that linelist (.cal)
        delta sigma / sigma against wavenumber plot (.eps)
        data used to make the plot (.plt)
    
    """
    print(docstring)


def calc_stat_unc(FWHM, del_w, SNR):
    return np.sqrt(FWHM * del_w) / SNR  # Same maths as = FWHM / (np.sqrt(FWHM / del_w) * SNR) = FWHM / (np.sqrt(Number of points across FWHM) * SNR)
    #return FWHM / (2 * SNR)
    
def get_header(header_lines_in):
    """Extracts the header lines at the start of aln and cln files.
    Returns them in a single string"""    
    
    header = ""
    
    for line in header_lines_in[:4]:
        header += line
        
    return header
    

def get_lines(get_lines_in):
    """Assigns each line in an aln or cln file to a Line object. 
    Returns a list of Line objects"""
    
    all_lines = []
    
    for line in get_lines_in[4:]:
        extracted_line = Line(line)  # create object of class Line
        all_lines.append(extracted_line)
         
    return all_lines
        
        
def remove_bad_lines(rm_bad_line_list, rm_bad_snr_disc_min, rm_bad_snr_disc_max):
    """Removes Line objects from line_list if they are not lsqfitted and/or
    have peak heights less than snr_discrim"""
    
    good_lines = []
    
    for line in rm_bad_line_list:  # only keep if lsqfitted and snr > discrinator
        if line.tags == "L" and float(line.peak) >= rm_bad_snr_disc_min and float(line.peak) <= rm_bad_snr_disc_max:
            good_lines.append(line)
    
    return good_lines
    
    
def get_common_lines(get_cmn_aln_lines, get_cmn_cln_lines, get_cmn_wavenum_discrim):
    """Returns list of tuples of (aln line, cln_line, delta sigma / sigma for that pair)
    that match lines in the cln line_list within a tolerance of +- wavenum_discrim"""
    
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
    

def get_corr_factor(get_corr_common_lines, std_dev_discrim, prev_calib_unc, aln_del_w, cln_del_w):
    """Discards lines outside of the std_dev_discrim limits and recalculates
    correction factor, std dev and std error for the remaining lines. Returns
    lines used in the new calcuation, the rejected lines, newly calculated
    correction factor, std dev and std error."""
     
    lines_used = []
    lines_rejected = [] 
    
    corr_factor = calc_corr_factor(get_corr_common_lines, prev_calib_unc, aln_del_w, cln_del_w)  # x used as dummy array as "weights" not needed here
    std_dev, std_error = calc_std_dev(get_corr_common_lines)
    limit = std_dev_discrim * std_dev  # set std_dev limit for delta sigma over sigma from corr_factor
    std_dev_limit = (corr_factor - limit, corr_factor + limit)
    
    for line in get_corr_common_lines: 
        delta_sig = line[2]  # delta sigma over sigma value for that alc/cln line pair
        if delta_sig <= std_dev_limit[1] and delta_sig >= std_dev_limit[0]:
            lines_used.append(line)
        else:
            lines_rejected.append(line)    
    
    corr_factor = calc_corr_factor(lines_used, prev_calib_unc, aln_del_w, cln_del_w)
    calib_unc = calc_calib_unc(lines_used, prev_calib_unc, aln_del_w, cln_del_w)
     
    return lines_used, lines_rejected, corr_factor, calib_unc, std_dev_limit


def calc_corr_factor(calc_corr_common_lines, prev_calib_unc, aln_del_w, cln_del_w):
    """Returns the weighted average of delta sigma over sigma for all good,
    common lines. This correction factor is weighted against SNR of the 
    lines"""
    
    weighted_factors = 0.
    total_weight = 0.
    
    for line in calc_corr_common_lines:
        ki = line[2]  # delta sigma over sigma
        aln_snr = float(line[0].peak)         
        aln_FWHM = float(line[0].width) / 1000  # to convert from mK to cm-1
        aln_stat_unc = calc_stat_unc(aln_FWHM, aln_del_w, aln_snr)  # Stat Unc. for aln line
        
        cln_snr = float(line[1].peak)
        cln_FWHM = float(line[1].width) / 1000
        cln_wavenum = float(line[1].wavenumber)
        cln_stat_unc = calc_stat_unc(cln_FWHM, cln_del_w, cln_snr)  # Stat Unc. for aln line
        cln_cal_unc = prev_calib_unc * cln_wavenum
        cln_tot_unc = np.sqrt(np.power(cln_stat_unc, 2) + np.power(cln_cal_unc, 2))  # prev_calib_unc is delta sigma / sigma therefore need to multiply by wavenumber
                
        wi = 1 / (np.power(aln_stat_unc, 2) + np.power(cln_tot_unc, 2))  # weighting for weighted average calc (inverse square of the aln and cln uncertainties added inquadrature)
        
        weighted_factors += (ki * wi)  # save total weighted factors for weighted avg calc
        total_weight += wi        
             
    return (weighted_factors / total_weight)


def calc_calib_unc(corr_fact_unc_common_lines, prev_calib_unc, aln_del_w, cln_del_w):
    """Calculates the calibration uncertainty. The calibration uncertainty for 
    the previous spectra is added to the current calibration uncertainty 
    (NOT in quadrature) to ensure that the calibration uncertainty increases 
    with distance from teh Ar II standards."""
    
    total_stat_unc = 0.
    snr_limit = 100. # Limit on SNR for strong lines to ensure sensible calibration uncertainty
        
    for line in corr_fact_unc_common_lines:
        aln_line = line[0]
        cln_line = line[1]        
        aln_SNR = float(aln_line.peak)
        aln_FWHM = float(aln_line.width) / 1000  # convert mK to cm-1                  
        cln_SNR = float(cln_line.peak)
        cln_FWHM = float(cln_line.width) / 1000  # convert mK to cm-1
               
        if aln_SNR > snr_limit:  
            aln_SNR = snr_limit
        if cln_SNR > snr_limit:
            cln_SNR = snr_limit

        aln_stat_unc = calc_stat_unc(aln_FWHM, aln_del_w, aln_SNR)  
        cln_stat_unc = calc_stat_unc(cln_FWHM, cln_del_w, cln_SNR)
               
        total_stat_unc += np.sqrt(np.power(aln_stat_unc, 2) + np.power(cln_stat_unc, 2))
                
    return prev_calib_unc + ((total_stat_unc) * np.power(10., -7))  #*10^-7 as this is the unc. on Keff which is in units of 10^-7 m
   
    
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


def plot_corr_factor(plot_lines_used, plot_lines_rejected, plot_corr_factor, plot_std_dev_limit, aln_filename, cln_filename, calib_unc):
    """Plots the lines used in (and lines rejected from) the correction factor 
    and calibration uncertainty calculations. The calculated correction factor 
    and the limits of the lines selection are also plotted."""
    
    lines_used_wn = []
    lines_used_del_sig = []
    lines_reject_wn = []
    lines_reject_del_sig = []  
    
    plot_data_filename = os.getcwd() + '/' + aln_filename.split('.')[0] + '.plt'
    plot_data = open(plot_data_filename, 'w')
    
    plot_data.write('Lines used in correction factor calculation:\n')
    plot_data.write('Wavenumber (cm-1) \t del sigma / sigma \t aln SNR \t cln SNR\n')
    
    for line in plot_lines_used:
        lines_used_wn.append(float(line[0].wavenumber))
        lines_used_del_sig.append(line[2])
        plot_data.write('%s \t %E \t %.1f \t %.1f\n' % (float(line[0].wavenumber), line[2], float(line[0].peak), float(line[1].peak)))
        
    plot_data.write('\n--------------------------------------------------\n')
    plot_data.write('Lines rejected from correction factor calcuation:\n')
    plot_data.write('Wavenumber (cm-1) \t del sigma / sigma\n')
    
    for line in plot_lines_rejected:
        lines_reject_wn.append(float(line[0].wavenumber))
        lines_reject_del_sig.append(line[2])
        plot_data.write('%s \t %E \t %.1f \t %.1f\n' % (float(line[0].wavenumber), line[2], float(line[0].peak), float(line[1].peak)))
    
    plot_data.write('\n--------------------------------------------------\n')
    plot_data.write('Correction factor = \t%E\n' % plot_corr_factor)
    plot_data.write('Upper s.d. limit  = \t%E\n' % plot_std_dev_limit[0])
    plot_data.write('Lower s.d. limit  = \t%E\n' % plot_std_dev_limit[1])
    plot_data.write('Calibration unc.  = \t %E' % calib_unc)
    plot_data.close()
     
    lines_used_del_sig = [x * np.power(10, 7) for x in lines_used_del_sig]  # convert everything into units of 10^-7
    lines_reject_del_sig = [x * np.power(10, 7) for x in lines_reject_del_sig]
    plot_corr_factor = plot_corr_factor * np.power(10, 7)
    plot_std_dev_limit = [x * np.power(10, 7) for x in plot_std_dev_limit]
    
    z = np.polyfit(lines_used_wn, lines_used_del_sig, 1)  # to calc and display linear trendline fot lines used.
    p = np.poly1d(z)        
   
    plt.plot(lines_used_wn, lines_used_del_sig, 'b+', label='Lines used (' + str(len(lines_used_wn)) + ')')
    plt.plot(lines_reject_wn, lines_reject_del_sig, 'r+', label='Lines rejected (' + str(len(lines_reject_wn)) + ')')
    plt.axhline(y=plot_corr_factor, color='b', linestyle='-')
    plt.axhline(y=(plot_std_dev_limit[0]), color='k', linestyle='--')
    plt.axhline(y=(plot_std_dev_limit[1]), color='k', linestyle='--')
    
    plt.ylabel(r'$\Delta \sigma / \sigma$ (10$^{-7}$)')
    plt.xlabel('Wavenumber (cm'r'$^{-1})$')
    plt.title(aln_filename + ' calibrated to ' + cln_filename)
    plt.legend()    
    plt.savefig(aln_filename.split('.')[0] + '.eps', format='eps', dpi=5000)
    
    plt.plot(lines_used_wn, p(lines_used_wn), 'b--')
    plt.show()



def write_cln_file(aln_lines, corr_factor, aln_header, aln_filename):
    """Outputs the wavenumber corrected xGremlin-format linelist file."""
    
    new_cln_filename = os.getcwd() + '/' + aln_filename.split('.')[0] + '.cln' 
    new_cln_file = open(new_cln_filename, 'w')
    
    aln_header = aln_header.split('0.0000000000000000')
    aln_header = aln_header[0] + str(corr_factor) + aln_header[1]
    new_cln_file.write(aln_header)
    
    for line in aln_lines[4:]:
        corrected_wavenum = '{0:.6f}'.format(float(line[8:20]) * (1 + corr_factor)) # correct wavenumber
        
        if float(corrected_wavenum) < 10000.0:  # to correct layout of cln files for 4 digit wavenumbers
            space_corr = ' '
        else:
            space_corr = ''
            
        new_cln_file.write(line[:8] + space_corr + str(corrected_wavenum) + line[20:])
    
    new_cln_file.close()
    

def write_cal_file(aln_raw_lines,  corr_factor, calib_unc, aln_filename, aln_del_w):
    """Ouputs a calibration file containing the lines written in the new
    calibrated linelist file aling with their line number and uncertainties.
    The total uncertainty is calculated as the addition of the statistical 
    unc. of the line and the calibration unc. of the spectrum in quadrature."""
    
    new_cal_filename = os.getcwd() + '/' + aln_filename.split('.')[0] + '.cal' 
    new_cal_file = open(new_cal_filename, 'w')
    
    new_cal_file.write('Correction factor = ' + '{0:.4E}'.format(corr_factor) + '\n')
    new_cal_file.write('Calibration uncertainty = ' + '{0:.4E}'.format(calib_unc) + '\n')
    new_cal_file.write('Line \tWavenumber \tStatisitcal Unc. \tTotal Unc. \n')
    
    for line in aln_raw_lines:
        
        corrected_wavenum = float(line.wavenumber) * (1 + corr_factor)
        
        if float(line.peak) > 200.:  # so as not to give an unreasonable estimate of uncertainty.
            line_snr = 200.
        else:
            line_snr = float(line.peak)
        
        FWHM = float(line.width) / 1000
        stat_unc = calc_stat_unc(FWHM, aln_del_w, line_snr)
        tot_unc = np.sqrt(np.power(stat_unc, 2) + np.power(calib_unc * corrected_wavenum, 2))
        stat_unc = '{0:.6f}'.format(stat_unc)
        tot_unc = '{0:.6f}'.format(tot_unc)      
        
        new_cal_file.write(line.line + '\t')
        new_cal_file.write('{0:.6f}'.format(corrected_wavenum) + '\t')
        new_cal_file.write(stat_unc + '\t')
        new_cal_file.write(tot_unc + '\n')
        
    new_cal_file.close()

        
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
        
class Decor:  # Here to make the terminal output look fancy!
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m' 
    

#### Main Program ####
def main():
    
    try:
        inp_file = open(sys.argv[1], 'r')    
    except IndexError:
        print_docstring()
        sys.exit()

#    inp_file = open('test.inp', 'r')  # For testing only 
    inp = inp_file.readlines()
    
    try: 
        snr_discrim_max = float(inp[4].split('=')[1])
        snr_discrim_min = float(inp[5].split('=')[1])
        wavenum_discrim = float(inp[6].split('=')[1])
        std_dev_discrim = float(inp[7].split('=')[1])
        prev_calib_unc = float(inp[8].split('=')[1])
        aln_del_w = float(inp[9].split('=')[1])
        cln_del_w = float(inp[10].split('=')[1])
        aln_filename = (inp[11].split('=')[1].strip())
        cln_filename = (inp[12].split('=')[1].strip()) 
    except IndexError:
        print('Input file is empty or corrupted. Please check and try again.')
        sys.exit()

    print('\nCalibrating %s against %s.\n' % (aln_filename, cln_filename))
    
    aln_in = open(aln_filename, 'r')
    cln_in = open(cln_filename, "r")
    
    aln_in_lines = aln_in.readlines()
    cln_in_lines = cln_in.readlines()

    aln_header = get_header(aln_in_lines)
    
    aln_raw_lines = get_lines(aln_in_lines)
    cln_raw_lines = get_lines(cln_in_lines)

    aln_good_lines = remove_bad_lines(aln_raw_lines, snr_discrim_min, snr_discrim_max)
    cln_good_lines = remove_bad_lines(cln_raw_lines, snr_discrim_min, snr_discrim_max)
    
    print('%s contains %i lsqfitted lines with a SNR between %.1f and %.1f.' % (aln_filename, len(aln_good_lines), snr_discrim_min, snr_discrim_max))
    print('%s contains %i lsqfitted lines with a SNR between %.1f and %.1f.\n' % (cln_filename, len(cln_good_lines), snr_discrim_min, snr_discrim_max))
    
    cmn_lines = get_common_lines(aln_good_lines, cln_good_lines, wavenum_discrim)
    
    print((Decor.UNDERLINE + Decor.BOLD + 'Lines matching within +- %.2f cm-1:'  + Decor.END) % wavenum_discrim)
    print('Line\taln Wavenumber\taln SNR \tcln SNR')
    
    for line in cmn_lines:
        print('%d\t%.4f \t%.1f    \t%.1f' % (float(line[0].line), float(line[0].wavenumber), float(line[0].peak), float(line[1].peak)))
    
    lines_used, lines_rejected, corr_factor, calib_unc, std_dev_limit = get_corr_factor(cmn_lines, std_dev_discrim, prev_calib_unc, aln_del_w, cln_del_w)
    
    print('\n\nRemoving lines from correction factor calculation based on distance from mean:')
    
    if len(lines_rejected) > 0:
        print((Decor.UNDERLINE + Decor.BOLD + '\nLines rejected from fit (> %s s.d. from mean):' + Decor.END) % std_dev_discrim)
    
        for line in lines_rejected:
            print('%d\t%.4f' % (float(line[0].line), float(line[0].wavenumber)))
    else:
        print('All lines within %s s.d. of mean.' % std_dev_discrim)
    
    print('\n\nCalculating final correction factor and calibration uncertainty:')
    print('\n------------------------------------------')
    print('Correction Factor       : %.4E' % corr_factor)
    print('Calibration Uncertainty : %.4E' % calib_unc)
    print('------------------------------------------\n')    
    
    plot_corr_factor(lines_used, lines_rejected, corr_factor, std_dev_limit, aln_filename.split('/')[-1], cln_filename.split('/')[-1], calib_unc)
    
    write_cln_file(aln_in_lines, corr_factor, aln_header, aln_filename.split('/')[-1])
    write_cal_file(aln_raw_lines, corr_factor, calib_unc, aln_filename.split('/')[-1], aln_del_w)
    
    aln_in.close()
    cln_in.close()   

if __name__ == '__main__':
    main()
