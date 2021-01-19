clc; clearvars; close all
% Script to display some MEG data on a dewar
% Currently requires brewer_maps to function and is only compatible with
% the CTF275 system. MEG data must be in an SPM D object

addpath('/Users/rtimms/Downloads/osl/brewer_maps')

%%% OPTIONS:
show_iso_lines=1;
N_isolines=1;
disp_anat=0;
tpts_of_interest=350:360; % what data to show from the D object. If a vector is provided then a movie will play
col_map='RdBu';




% Load in some co-registered MEG data
% D_coreg=spm_eeg_load( '/Volumes/TASER/Notts_motor/03677_processed/03677341_Ellie_20170615_04/edff03677341_Ellie_20170615_04.mat')
D_coreg=spm_eeg_load('/Volumes/TASER/Notts_motor/11251_processed/11251_TAP_08/edff11251_TAP_08.mat');



show_data_on_dewar(D_coreg,D_coreg(D_coreg.indchantype('MEEG','GOOD'),1,1),show_iso_lines,disp_anat,N_isolines,tpts_of_interest,col_map)

