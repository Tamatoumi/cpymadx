#-----------------------------------------------------------------------------------------
# Tatiana Da Silva - CEA
#-----------------------------------------------------------------------------------------

from cpymad.madx import Madx
import matplotlib.pyplot as plt
from scipy.stats import norm
import def_functions as df
import numpy as np
import pandas as pd
import sys
import os


madx = Madx()

madx.option(echo=True)

os.system("mkdir -p temp")
os.system("mkdir -p results")

madx.set(format="12d")
madx.set(format="20.12g")
madx.set(format="-18s")

madx.option(rbarc=False)

#----- SEQUENCE DEFINITION -----
madx.call(file='/home/td271008/work/cpymadx/FCCee_heb_modett.ele')
madx.call(file='/home/td271008/work/cpymadx/FCCee_heb_modett.seq')
madx.call(file='/home/td271008/work/cpymadx/FCCee_heb_modett.str')

#----- MACRO DEFINITION -----
madx.call(file='/home/td271008/work/madx/FCCee_heb_misc_macros.madx')
madx.call(file='/home/td271008/work/errors/FCCee_errors_macros.madx')

beam_mode=3
injection_mode=1
with_radiate=0
madx.exec('load_beam()')
madx.exec('tws_select()')
tol_tar=1e-12
switch_rf_on=1
vrf400=60  #RF Voltage at 20 GeV
vrf800=0

#----- REFERENCE OPTICS -----
madx.command.beam(sequence='fcc_heb', particle='electron')
madx.use(sequence='fcc_heb')
madx.seqedit(sequence='fcc_heb')
madx.flatten()
madx.cycle(start='ip1')
madx.endedit()
madx.use(sequence='fcc_heb')
madx.savebeta(label='BETAIP1',place='IP1')
madx.twiss(sequence='fcc_heb', file='FCCee_heb_modett.tfs')
#fcc_heb=madx.sequence.fcc_heb
#twiss11=fcc_heb.twiss
twiss1=madx.twiss(sequence='fcc_heb')
#print(twiss1.x)
#sys.exit()
'''
MQ_A=df.get_elem(madx.elements,'mq','a')
CORR_A=df.get_elem(madx.elements,'corr','a')
BPM_A=df.get_elem(madx.elements,'bpm','a')
'''
#print(twiss1.betx)
#sys.exit()

summ=madx.table.summ
qx_tar=summ.q1
qy_tar=summ.q2
dqx_tar=0
dqy_tar=0
          
#----- ERRORS -----
madx.call(file='/home/td271008/work/optics_RF_9090/select_errors.madx')

#----- SAVE ERRORS ----
madx.option(echo=True)
madx.select(flag='error')
madx.esave(file='FCCee_heb_errors_corr.out')

#----- ANALYTICAL CORRECTORS CALCULATION -----
file1="FCCee_heb_modett.tfs"
file2="FCCee_heb_errors_corr.out"

df.anal_corr_calc(file1,file2)

#----- CORRECT TUNES & CHROMA -----

madx.call(file='/home/td271008/work/cpymadx/analytical_corr.str')
#madx.call(file='/home/td271008/work/madx/FCCee_heb-quads.str') #do not use with optic modett

'''
#tune match
df.tune_match('fcc_heb',1.0E-7)

#chroma match
df.chroma_match('fcc_heb',1.0E-7)
'''

madx.exec('SEXTUOFF')

madx.twiss(sequence='fcc_heb', file='FCCee_heb_modett_err.tfs')

#twiss2=madx.twiss(sequence='fcc_heb')
madx.twiss(beta0='BETAIP1')

#----- CORRECTION -----

madx.command.usemonitor(status='on')
#madx.command.usekick(status='on')
madx.command.usekick(status='off')
madx.command.usekick(status='on', range='S.DIS_LEFT.A1/E.DIS_RIGHT.A1')
madx.command.usekick(status='on', range='S.DIS_LEFT.A2/E.DIS_RIGHT.A2')
madx.command.usekick(status='on', range='S.DIS_LEFT.A3/E.DIS_RIGHT.A3')
madx.command.usekick(status='on', range='S.DIS_LEFT.A4/E.DIS_RIGHT.A4')
madx.command.usekick(status='on', range='S.DIS_LEFT.A5/E.DIS_RIGHT.A5')
madx.command.usekick(status='on', range='S.DIS_LEFT.A6/E.DIS_RIGHT.A6')
madx.command.usekick(status='on', range='S.DIS_LEFT.A7/E.DIS_RIGHT.A7')
madx.command.usekick(status='on', range='S.DIS_LEFT.A8/E.DIS_RIGHT.A8')

#madx.command.correct(flag='line',mode='svd',monerror=1,error=1.0e-7,plane='x',clist='results/cx_fccee_heb.tab')
#madx.command.correct(flag='line',mode='svd',monerror=1,error=1.0e-7,plane='y',clist='results/cy_fccee_heb.tab')

madx.command.correct(flag='line',mode='svd',monerror=1,error=1.0e-7,plane='x')
px_corr_x=list(madx.table.corr['px.correction'])
py_corr_x=list(madx.table.corr['py.correction'])

madx.command.correct(flag='line',mode='svd',monerror=1,error=1.0e-7,plane='y')
px_corr_y=list(madx.table.corr['px.correction'])
py_corr_y=list(madx.table.corr['py.correction'])


#madx.twiss(sequence='fcc_heb', file='FCCee_heb_modett_orbcor.tfs')
#madx.twiss(BETA0='BETAIP1', sequence='fcc_heb', file='FCCee_heb_4IPs_ttRF_9090_qx565_qy595_mqonly_orbcor.tfs')
twiss3=madx.twiss(beta0='BETAIP1', sequence='fcc_heb')
#print(twiss3.betx)



'''
madx.command.correct(flag='line',mode='svd',monerror=1,error=1.0e-7,plane='x',clist='results/cx_fccee_heb_it2.tab')
madx.command.correct(flag='line',mode='svd',monerror=1,error=1.0e-7,plane='y',clist='results/cy_fccee_heb_it2.tab')

madx.twiss(BETA0='BETAIP1', sequence='fcc_heb', file='FCCee_heb_4IPs_ttRF_9090_qx565_qy595_mqonly_orbcor_it2.tfs')
'''

#----- PLOT -----


''''
print('\n')
print('\n')
#print(twiss3.betx-twiss1.betx)
print('\n')
'''

df.plot_graph(twiss1,twiss3)
sys.exit()


plt.rcParams.update({'font.size':13})
plt.rc('xtick',labelsize=11)
plt.rc('ytick',labelsize=11)

fnameref="FCCee_heb_4IPs_ttRF_9090_qx565_qy595.tfs"
head_opt=pd.read_csv(fnameref, header=50, sep='\s+', nrows=0).columns[1:]
optics_ref=pd.read_csv(fnameref, skiprows=52, sep='\s+', names=head_opt)
print(optics_ref)

fnamerr='FCCee_heb_4IPs_ttRF_9090_qx565_qy595_mqonly_orbcor.tfs'
head_opt=pd.read_csv(fnamerr, header=50, sep='\s+', nrows=0).columns[1:]
optics_err=pd.read_csv(fnamerr, skiprows=52, sep='\s+', names=head_opt)
optics_err = optics_err.reset_index(drop=True)
print(optics_err)

print('\n')
print('\n')
print(optics_err["BETX"]-optics_ref["BETX"])
print('\n')
print('\n')

# Plot orbits in m
fig, ax=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax[0].plot(optics_err["S"]/1000., optics_err["X"], "-r")
ax[0].set_ylabel("x [m]")
ax[1].plot(optics_err["S"]/1000., optics_err["Y"], "-b")
ax[1].set_xlabel("longitundinal position [km]")
ax[1].set_ylabel("y [m]")
ax[1].set_xlim(0,93)
fig, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
plt.title("Orbits in m")

# Plot beta-beat in %
fig1, ax1=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax1[0].plot(optics_err["S"]/1000., 100.*((optics_err["BETX"]-optics_ref["BETX"])/optics_ref["BETX"]), "-r")
ax1[0].set_ylabel(r"$\Delta\beta_{x}/\beta_{x,ref}$ [%]")
ax1[1].plot(optics_err["S"]/1000., 100.*((optics_err["BETY"]-optics_ref["BETY"])/optics_ref["BETY"]), "-b")
ax1[1].set_xlabel("longitundinal position [km]")
ax1[1].set_ylabel(r"$\Delta\beta_{y}/\beta_{y,ref}$ [%]")
ax1[1].set_xlim(0,93)
fig1, plt.subplots_adjust(left=.13, right=.97, top=.94, bottom=.11)
plt.title("Beta-beat in %")

# Plot normalized dispersion mm^1/2
fig2, ax2=plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
ax2[0].plot(optics_err["S"]/1000., (optics_err["DX"]-optics_ref["DX"])/np.sqrt(optics_ref["BETX"]), "-r")
ax2[0].set_ylabel(r"$\Delta D_{x}/\sqrt{\beta_{x,ref}}$ [m$^{1/2}$]")
ax2[1].plot(optics_err["S"]/1000., (optics_err["DY"]-optics_ref["DY"])/np.sqrt(optics_ref["BETY"]), "-b")
ax2[1].set_xlabel("longitundinal position [km]")
ax2[1].set_ylabel(r"$\Delta D_{y}/\sqrt{\beta_{y,ref}}$ [m$^{1/2}$]")
ax2[1].set_xlim(0,93)
fig2, plt.subplots_adjust(left=.16, right=.97, top=.94, bottom=.11)
plt.title("Normalized dispersion [mm$^{1/2}$]")


plt.show()
