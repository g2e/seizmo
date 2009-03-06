Installation Procedure


The installation of TTBOX is relatively simple. All you have to do is to create some
directories and tell MatLab, where TTBOX is. 



1 if it does not already exist, create a directory matlab within your home directory. 

2 Copy the TTBOX directory into the matlab directory 

3 if it does not already exist, create a file startup.m within the matlab directory 

4 Insert the following lines into startup.m: 

  basepath='This should be the path to your home directory!';
  addpath([basepath filesep 'matlab' filesep 'ttbox'],1);
  addpath([basepath filesep 'matlab' filesep 'ttbox' filesep 'support'],1);
      
5 execute startup within the MatLab command window. 



MatLab is now ready to execute TTBOX routines. 
The startup-script will be executed automatically upon start of MatLab. 

In order to read the online documentation, you have to browse to the doc directory
and open the ttbox/doc/ttbox.html file with your web browser (it is recommended to
set a bookmark). 
It is not important where the TTBOX routines and the support directory reside.
What is important, is that an addpath-call during startup tells MatLab where it 
actually is. Please refer to the MatLab documentation to find out about startup 
behavior of MatLab and its search path list. 

