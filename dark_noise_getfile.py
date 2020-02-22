###################################
#SCUBA-2 Dark Noise Comparisons
###################################
'''
INPUT: Path to log files on EAO computer
       Output Directory
OUTPUT: Effective NEP vs Time plots for 450/850 um 
        Histogram plots of Effective NEP split by sub-array for 450/850 um
NOTES: Needs to be run on EAO computer.

Last Updated: Feb 22, 2020
Written by: A.J. Tetarenko, based on D. Bintley's C-shell script
'''
###################################

import numpy as np
import glob
import os
from astropy.time import Time
from astropy.io import ascii
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.ticker import AutoMinorLocator
import matplotlib.dates as mdates

def grab_files(path_to_logs,file_out):
	''' Digs through the log files and grabs the Effective NEPs and dates/times
	for each sub-array, then writes them to separate output files for 450/850 um'''
	nep_eff8=[]
	times8=[]
	subarr8=[]
	nep_eff4=[]
	times4=[]
	subarr4=[]
	skip4=[]
	skip8=[]
	dates=glob.glob(path_to_logs+'*')
	print(dates[-1])
	for i in range(0,len(dates)):
		path_dir=dates[i]
		day=path_dir.split('/')[-1]
		if os.path.isfile(path_dir+'/850-summit/log.nep'):
			try:
				log_ncol=ascii.read(path_dir+'/850-summit/log.nep')
				if len(log_ncol.columns)==11:
					log_nep8=ascii.read(path_dir+'/850-summit/log.nep',names=('UT','HST','Obs','Subscan','Subarray','Median','Mean','Nbol','NEP_eff','NEP_wt','Shutter'))
					#log_bnoi8=ascii.read(path_dir+'/850-summit/log.bolonoise',names=('UT','HST','Obs','Subscan','Subarray','Median','Mean','Nbol','NEP_eff','NEP_wt','Shutter'))
				elif np.all(log_ncol['col12']<1) and np.all(log_ncol['col12']>0):
					log_nep8=ascii.read(path_dir+'/850-summit/log.nep',names=('UT','HST','Obs','Subscan','Subarray','Median','Mean','Nbol','NEP_eff','NEP_wt','Shutter','Sigma'))
					#log_bnoi8=ascii.read(path_dir+'/850-summit/log.bolonoise',names=('UT','HST','Obs','Subscan','Subarray','Median','Mean','Nbol','NEP_eff','NEP_wt','Shutter','Sigma'))
				else:
					log_nep8=ascii.read(path_dir+'/850-summit/log.nep',names=('UT','HST','Obs','Subscan','Subarray','Median','Mean','Mode','Nbol','NEP_eff','NEP_wt','Shutter'))
					#log_bnoi8=ascii.read(path_dir+'/850-summit/log.bolonoise',names=('UT','HST','Obs','Subscan','Subarray','Median','Mean','Nbol','NEP_eff','NEP_wt','Shutter'))
				nep_eff8.extend(log_nep8['NEP_eff'][log_nep8['Shutter']==0.0])
				times8.extend(log_nep8['HST'][log_nep8['Shutter']==0.0])
				subarr8.extend(log_nep8['Subarray'][log_nep8['Shutter']==0.0])
			except Exception:
				foo=2
				skip8.append(day)
		if os.path.isfile(path_dir+'/450-summit/log.nep'):
			try:
				log_ncol=ascii.read(path_dir+'/450-summit/log.nep')
				if len(log_ncol.columns)==11:
					log_nep4=ascii.read(path_dir+'/450-summit/log.nep',names=('UT','HST','Obs','Subscan','Subarray','Median','Mean','Nbol','NEP_eff','NEP_wt','Shutter'))
				elif np.all(log_ncol['col12']<1) and np.all(log_ncol['col12']>0):
					log_nep4=ascii.read(path_dir+'/450-summit/log.nep',names=('UT','HST','Obs','Subscan','Subarray','Median','Mean','Nbol','NEP_eff','NEP_wt','Shutter','Sigma'))
					#log_bnoi4=ascii.read(path_dir+'/450-summit/log.bolonoise',names=('UT','HST','Obs','Subscan','Subarray','Median','Mean','Nbol','NEP_eff','NEP_wt','Shutter','Sigma'))
				else:
					log_nep4=ascii.read(path_dir+'/450-summit/log.nep',names=('UT','HST','Obs','Subscan','Subarray','Median','Mean','Mode','Nbol','NEP_eff','NEP_wt','Shutter'))
					#log_bnoi4=ascii.read(path_dir+'/450-summit/log.bolonoise',names=('UT','HST','Obs','Subscan','Subarray','Median','Mean','Nbol','NEP_eff','NEP_wt','Shutter'))
				nep_eff4.extend(log_nep4['NEP_eff'][log_nep4['Shutter']==0.0])
				times4.extend(log_nep4['HST'][log_nep4['Shutter']==0.0])
				subarr4.extend(log_nep4['Subarray'][log_nep4['Shutter']==0.0])
			except Exception:
				foo=2
				skip4.append(day)
	fileo=open(file_out+'s8.txt','w')
	for i in range(0,len(times8)):
		if '-' in times8[i] and nep_eff8[i]!=0.0:
			fileo.write('{0} {1} {2} {3}\n'.format(Time(times8[i],format='isot').mjd,times8[i],nep_eff8[i],subarr8[i]))
	fileo.close()
	fileo=open(file_out+'s4.txt','w')
	for i in range(0,len(times4)):
		if '-' in times4[i] and nep_eff4[i]!=0.0:
			fileo.write('{0} {1} {2} {3}\n'.format(Time(times4[i],format='isot').mjd,times4[i],nep_eff4[i],subarr4[i]))
	fileo.close()

###################################
#USER INPUT
###################################
path_to_logs='/jcmtdata/raw/pipeline-log/'
file_out='/export/data2/atetarenko/dark_noise/files/'

#event dates to add to time series plots
bad_dates={'Mem. Off': '2017-12-06', 'Mem. On': '2018-04-30',\
'New Filters': '2016-11-19', 'SMU Fix': '2018-07-27'}

#start date for Effective NEP histograms
date_after='2018-01-01'
###################################

#grab the info needed from JCMT logs
grab_files(path_to_logs,file_out)

#read in files created above
s8=ascii.read(file_out+'s8.txt',names=('MJD','ISOT','NEP','Subarr'))
s4=ascii.read(file_out+'s4.txt',names=('MJD','ISOT','NEP','Subarr'))

#Effective NEP vs time for 850um
plt.rcdefaults()
font={'family':'serif','size':'14'}
rc('font',**font)
mpl.rcParams['xtick.direction'] ='in'
mpl.rcParams['ytick.direction'] ='in'
fig=plt.figure()
ax=plt.subplot(111)
#ax2=ax.twiny()
ax.plot(Time(s8['MJD'][s8['Subarr']=='s8a'],format='mjd').datetime, s8['NEP'][s8['Subarr']=='s8a']*1e18,label='s8a',color='b',marker='o',ls='',ms=2)#,markevery=10)
ax.plot(Time(s8['MJD'][s8['Subarr']=='s8b'],format='mjd').datetime, s8['NEP'][s8['Subarr']=='s8b']*1e18,label='s8b',color='orange',marker='o',ls='',ms=2)#,markevery=10)
ax.plot(Time(s8['MJD'][s8['Subarr']=='s8c'],format='mjd').datetime, s8['NEP'][s8['Subarr']=='s8c']*1e18,label='s8c',color='g',marker='o',ls='',ms=2)#,markevery=10)
ax.plot(Time(s8['MJD'][s8['Subarr']=='s8d'],format='mjd').datetime, s8['NEP'][s8['Subarr']=='s8d']*1e18,label='s8d',color='r',marker='o',ls='',ms=2)#,markevery=10)
for item in bad_dates.keys():
	ax.axvline(x=Time(bad_dates[item],format='iso').datetime, color='k',ls=':')
	ax.text(Time(bad_dates[item],format='iso').datetime,1500, item,rotation=90,fontsize=10)
ax.set_xlabel('Time')
ax.set_ylabel('Effective NEP ($\\times 10^{-18}$)')
ax.legend(loc='upper left', ncol=2,facecolor='w',framealpha=1)
ax.set_ylim(1,6000)
#ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.set_yscale('log')
years=mdates.YearLocator(1)
months=mdates.MonthLocator(7)
ax.xaxis.set_major_locator(years)
ax.xaxis.set_minor_locator(months)
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
plt.setp(ax.get_xticklabels(),rotation=45,horizontalalignment='right')
ax.tick_params(axis='x', which='major',labelsize=15,length=7, width=1.5,top='on',bottom='on',pad=7)
ax.tick_params(axis='x', which='minor',labelsize=15,length=5, width=1.,top='on',bottom='on',pad=7)
ax.tick_params(axis='y', which='major',labelsize=15,length=7, width=1.5,left='on',right='on',pad=7)
ax.tick_params(axis='y', which='minor',labelsize=15,length=5, width=1.,left='on',right='on',pad=7)
plt.savefig(file_out+'s850_NEPvstime.pdf',bbox_inches='tight')
plt.show()

#Effective NEP vs time for 450um
fig=plt.figure()
ax=plt.subplot(111)
#ax2=ax.twiny()
ax.plot(Time(s4['MJD'][s4['Subarr']=='s4a'],format='mjd').datetime, s4['NEP'][s4['Subarr']=='s4a']*1e18,label='s4a',color='b',marker='o',ls='',ms=2)#,markevery=10)
ax.plot(Time(s4['MJD'][s4['Subarr']=='s4b'],format='mjd').datetime, s4['NEP'][s4['Subarr']=='s4b']*1e18,label='s4b',color='orange',marker='o',ls='',ms=2)#,markevery=10)
ax.plot(Time(s4['MJD'][s4['Subarr']=='s4c'],format='mjd').datetime, s4['NEP'][s4['Subarr']=='s4c']*1e18,label='s4c',color='g',marker='o',ls='',ms=2)#,markevery=10)
ax.plot(Time(s4['MJD'][s4['Subarr']=='s4d'],format='mjd').datetime, s4['NEP'][s4['Subarr']=='s4d']*1e18,label='s4d',color='r',marker='o',ls='',ms=2)#,markevery=10)
ax.set_xlabel('Time')
ax.set_ylabel('Effective NEP ($\\times 10^{-18}$)')
ax.legend(loc='upper left', ncol=2,facecolor='w',framealpha=1)
for item in bad_dates.keys():
	ax.axvline(x=Time(bad_dates[item],format='iso').datetime, color='k',ls=':')
	ax.text(Time(bad_dates[item],format='iso').datetime,8000, item,rotation=90,fontsize=10)
#ax.set_ylim(3,200)
#ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.set_yscale('log')
years=mdates.YearLocator(1)
months=mdates.MonthLocator(7)
ax.set_ylim(1,50000)
ax.xaxis.set_major_locator(years)
ax.xaxis.set_minor_locator(months)
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
plt.setp(ax.get_xticklabels(),rotation=45,horizontalalignment='right')
ax.tick_params(axis='x', which='major',labelsize=15,length=7, width=1.5,top='on',bottom='on',pad=7)
ax.tick_params(axis='x', which='minor',labelsize=15,length=5, width=1.,top='on',bottom='on',pad=7)
ax.tick_params(axis='y', which='major',labelsize=15,length=7, width=1.5,left='on',right='on',pad=7)
ax.tick_params(axis='y', which='minor',labelsize=15,length=5, width=1.,left='on',right='on',pad=7)
plt.savefig(file_out+'s450_NEPvstime.pdf',bbox_inches='tight')
plt.show()



#Effective NEP Histograms for 850um
s8a=s8[s8['Subarr']=='s8a']
s8b=s8[s8['Subarr']=='s8b']
s8c=s8[s8['Subarr']=='s8c']
s8d=s8[s8['Subarr']=='s8d']
fig=plt.figure()
ax=plt.subplot(221)
ax1=plt.subplot(222)
ax2=plt.subplot(223)
ax3=plt.subplot(224)
ax.hist(np.array(s8a['NEP'][s8a['MJD']>Time(date_after,format='iso').mjd])*1e18,bins=np.logspace(np.log10(1),np.log10(200),50),label='s8a',alpha=0.7,color='b')
ax1.hist(np.array(s8b['NEP'][s8b['MJD']>Time(date_after,format='iso').mjd])*1e18,bins=np.logspace(np.log10(1),np.log10(200),50),label='s8b',alpha=0.7,color='orange')
ax2.hist(np.array(s8c['NEP'][s8c['MJD']>Time(date_after,format='iso').mjd])*1e18,bins=np.logspace(np.log10(1),np.log10(200),50),label='s8c',alpha=0.7,color='g')
ax3.hist(np.array(s8d['NEP'][s8d['MJD']>Time(date_after,format='iso').mjd])*1e18,bins=np.logspace(np.log10(1),np.log10(200),50),label='s8d',alpha=0.7,color='r')
ax.set_yscale('log')
ax.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xscale('log')
ax3.set_yscale('log')
ax3.set_xscale('log')
ax.tick_params(axis='x', which='major',labelsize=15,length=7, width=1.5,top='on',bottom='on',pad=7)
ax.tick_params(axis='x', which='minor',labelsize=15,length=5, width=1.,top='on',bottom='on',pad=7)
ax.tick_params(axis='y', which='major',labelsize=15,length=7, width=1.5,left='on',right='on',pad=7)
ax.tick_params(axis='y', which='minor',labelsize=15,length=5, width=1.,left='on',right='on',pad=7)
ax1.tick_params(axis='x', which='major',labelsize=15,length=7, width=1.5,top='on',bottom='on',pad=7)
ax1.tick_params(axis='x', which='minor',labelsize=15,length=5, width=1.,top='on',bottom='on',pad=7)
ax1.tick_params(axis='y', which='major',labelsize=15,length=7, width=1.5,left='on',right='on',pad=7)
ax1.tick_params(axis='y', which='minor',labelsize=15,length=5, width=1.,left='on',right='on',pad=7)
ax2.tick_params(axis='x', which='major',labelsize=15,length=7, width=1.5,top='on',bottom='on',pad=7)
ax2.tick_params(axis='x', which='minor',labelsize=15,length=5, width=1.,top='on',bottom='on',pad=7)
ax2.tick_params(axis='y', which='major',labelsize=15,length=7, width=1.5,left='on',right='on',pad=7)
ax2.tick_params(axis='y', which='minor',labelsize=15,length=5, width=1.,left='on',right='on',pad=7)
ax3.tick_params(axis='x', which='major',labelsize=15,length=7, width=1.5,top='on',bottom='on',pad=7)
ax3.tick_params(axis='x', which='minor',labelsize=15,length=5, width=1.,top='on',bottom='on',pad=7)
ax3.tick_params(axis='y', which='major',labelsize=15,length=7, width=1.5,left='on',right='on',pad=7)
ax3.tick_params(axis='y', which='minor',labelsize=15,length=5, width=1.,left='on',right='on',pad=7)
ax2.set_xlabel('Effective NEP ($\\times 10^{-18}$)')
ax2.xaxis.set_label_coords(1.1,-0.2)
ax2.set_ylabel('Number of Occurences')
ax2.yaxis.set_label_coords(-0.25,1.1)
ax.legend(loc='upper right', ncol=1)
ax1.legend(loc='upper right', ncol=1)
ax2.legend(loc='upper right', ncol=1)
ax3.legend(loc='upper right', ncol=1)
ax.set_title(date_after+' - Present',x=1.1,y=1)
ax.set_xticklabels([])
ax1.set_xticklabels([])
ax1.set_yticklabels([])
ax3.set_yticklabels([])
plt.subplots_adjust(hspace=0.1,wspace=0.1)
plt.savefig(file_out+'s850_NEP_hists.pdf',bbox_inches='tight')
plt.show()

#Effective NEP Histograms for 450um
s4a=s4[s4['Subarr']=='s4a']
s4b=s4[s4['Subarr']=='s4b']
s4c=s4[s4['Subarr']=='s4c']
s4d=s4[s4['Subarr']=='s4d']
fig=plt.figure()
ax=plt.subplot(221)
ax1=plt.subplot(222)
ax2=plt.subplot(223)
ax3=plt.subplot(224)
ax.hist(np.array(s4a['NEP'][s4a['MJD']>Time(date_after,format='iso').mjd])*1e18,bins=np.logspace(np.log10(1),np.log10(200),50),label='s4a',alpha=0.7,color='b')
ax1.hist(np.array(s4b['NEP'][s4b['MJD']>Time(date_after,format='iso').mjd])*1e18,bins=np.logspace(np.log10(1),np.log10(200),50),label='s4b',alpha=0.7,color='orange')
ax2.hist(np.array(s4c['NEP'][s4c['MJD']>Time(date_after,format='iso').mjd])*1e18,bins=np.logspace(np.log10(1),np.log10(200),50),label='s4c',alpha=0.7,color='g')
ax3.hist(np.array(s4d['NEP'][s4d['MJD']>Time(date_after,format='iso').mjd])*1e18,bins=np.logspace(np.log10(1),np.log10(200),50),label='s4d',alpha=0.7,color='r')
ax.set_yscale('log')
ax.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xscale('log')
ax3.set_yscale('log')
ax3.set_xscale('log')
ax.tick_params(axis='x', which='major',labelsize=15,length=7, width=1.5,top='on',bottom='on',pad=7)
ax.tick_params(axis='x', which='minor',labelsize=15,length=5, width=1.,top='on',bottom='on',pad=7)
ax.tick_params(axis='y', which='major',labelsize=15,length=7, width=1.5,left='on',right='on',pad=7)
ax.tick_params(axis='y', which='minor',labelsize=15,length=5, width=1.,left='on',right='on',pad=7)
ax1.tick_params(axis='x', which='major',labelsize=15,length=7, width=1.5,top='on',bottom='on',pad=7)
ax1.tick_params(axis='x', which='minor',labelsize=15,length=5, width=1.,top='on',bottom='on',pad=7)
ax1.tick_params(axis='y', which='major',labelsize=15,length=7, width=1.5,left='on',right='on',pad=7)
ax1.tick_params(axis='y', which='minor',labelsize=15,length=5, width=1.,left='on',right='on',pad=7)
ax2.tick_params(axis='x', which='major',labelsize=15,length=7, width=1.5,top='on',bottom='on',pad=7)
ax2.tick_params(axis='x', which='minor',labelsize=15,length=5, width=1.,top='on',bottom='on',pad=7)
ax2.tick_params(axis='y', which='major',labelsize=15,length=7, width=1.5,left='on',right='on',pad=7)
ax2.tick_params(axis='y', which='minor',labelsize=15,length=5, width=1.,left='on',right='on',pad=7)
ax3.tick_params(axis='x', which='major',labelsize=15,length=7, width=1.5,top='on',bottom='on',pad=7)
ax3.tick_params(axis='x', which='minor',labelsize=15,length=5, width=1.,top='on',bottom='on',pad=7)
ax3.tick_params(axis='y', which='major',labelsize=15,length=7, width=1.5,left='on',right='on',pad=7)
ax3.tick_params(axis='y', which='minor',labelsize=15,length=5, width=1.,left='on',right='on',pad=7)
ax2.set_xlabel('Effective NEP ($\\times 10^{-18}$)')
ax2.xaxis.set_label_coords(1.1,-0.2)
ax2.set_ylabel('Number of Occurences')
ax2.yaxis.set_label_coords(-0.25,1.1)
ax.legend(loc='upper right', ncol=1)
ax1.legend(loc='upper right', ncol=1)
ax2.legend(loc='upper right', ncol=1)
ax3.legend(loc='upper right', ncol=1)
ax.set_title(date_after+' - Present',x=1.1,y=1)
ax.set_xticklabels([])
ax1.set_xticklabels([])
ax1.set_yticklabels([])
ax3.set_yticklabels([])
plt.subplots_adjust(hspace=0.1,wspace=0.1)
plt.savefig(file_out+'s450_NEP_hists.pdf',bbox_inches='tight')
plt.show()

