# SPDC TYPE II bLOB
###################################
import math
import numpy as np
import matplotlib.pyplot as plt
import os, time, glob
#***********PARAMETERS****************
def SPDC(wavep, angle, L, distz, wp):
    wavep = wavep/1000.0  # wavelength of pump field in 405nm ### angle = 28.9     # Angle in degree
    thetap = np.radians(angle)  # thetap in radians45.64, 42.58 40.5 42.148
    L= L*1000.00      # Crystal thickness in 2mm
    distz= distz*10000.0      # distz in cm 100cm  , wp in  um
    ########################################
    #*******Dispersion relation o ray***************
    def noo(n):
        global wavep
        n2 = 2.7405 + ((0.0184)/((n**2)-0.0179)) - (0.0155 * n**2)
        rio = round(math.sqrt(n2),4)
        return rio
    #********Dispersion relation e ray***************
    def neo(n):
        global wavep
        n2 = 2.3730 + ((0.0128) / ((n ** 2) - 0.0156)) - (0.0044 * n ** 2)
        rio = round(math.sqrt(n2), 4)
        return rio
    #############################################
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #####Two photon wavefunction##############
    def two_photon_wavefunction(xs,ys,xi,yi,wp):
        ss= (np.square(xs)+ np.square(ys)) + (np.square(xi)+ np.square(yi))+ (2*np.sqrt( (np.square(xs)+ np.square(ys))*(np.square(xi)+ np.square(yi)))*np.cos(np.arctan2(ys,xs)-np.arctan2(yi,xi)))
        kp = (2 * (np.pi)) / wavep
        kpz = (kp * etap) + (alphap * (xs + xi)) - ((1 / (2 * kp * etap)) * (np.square(betap * (xs + xi)) + np.square(gammap * (ys + yi))))
        ksz = (kp * nobar / 2) - ((1 / (kp * nobar)) * (np.square(xs) + np.square(ys)))
        kiz = (kp * nobar / 2) - ((1 / (kp * nobar)) * (np.square(xi) + np.square(yi)))
        tmp = ((ksz + kiz - kpz) * (L / 2))
    #-------------------------------------------------------------------------------------
        ksze=((kp/2.0) * etas) + (alphas * (xs)) - ((1 / ( kp * etas)) * (np.square(betas * xs ) + np.square(gammas * ys)))
        tmpeo = ((ksze + kiz - kpz) * (L / 2))
    #----------------------------------------------------------------------------------------------
        kize=((kp/2.0) * etai) + (alphai * (xi)) - ((1 / ( kp * etai)) * (np.square(betai * xi) + np.square(gammai * yi)))
        tmpoe = ((ksz + kize - kpz) * (L / 2))
    #--------------------------------------------------------------------------------------------
        pumpfield =np.exp(-(wp**2 + ((2j)*(distz/(kp*etap))))*(ss/4.0))
        fullfunction =pumpfield*((np.sinc(tmpeo/(np.pi))*np.exp(-1.0j *tmpeo))+(np.sinc(tmpoe/(np.pi))*np.exp(-1.0j *tmpoe))+(np.sinc(tmp/(np.pi))*np.exp(-1.0j *tmp)))
        return fullfunction
    #-----------------------------------------------------
    #***************************************************
    if angle <= 41.0:
        rho = 0.6
    elif angle > 41 and angle < 41.5 :
        rho = 0.8
    elif angle >= 41.5 and angle < 41.8 :
        rho= 0.9
    elif angle > 41.8 and angle < 45 :
        rho= 1.0
    elif angle >= 45.0 and angle < 50.0 :
        rho= 1.5
    elif angle >= 50.0 and angle < 65.0 :
        rho = 2.0
    elif angle >=65 :
        rho = 2.8
    ###print (rho)
    ###rho=2.8
    grdpnt=200
    dx= (2*rho/grdpnt)
    xs,ys = np.meshgrid(np.linspace(-rho, rho,grdpnt), np.linspace(-rho,rho,grdpnt))
    # alpha ,beta ,gamma, eta ***********************
    no = noo(wavep)
    ne = neo(wavep)
    den = np.square(no * np.sin(thetap)) + np.square(ne * np.cos(thetap))
    alphap = ((np.square(no) - np.square(ne)) * (np.sin(thetap)) * (np.cos(thetap))) / den
    betap = (no * ne) / den
    gammap = no / math.sqrt(den)
    etap = ne * gammap
    nobar = noo(2 * wavep)
    #------------------8888888888888------------------------
    nos = noo(2*wavep)
    nes = neo(2*wavep)
    thetas= thetap - np.arctan2(ys,xs)
    dens = np.square(nos * np.sin(thetap)) + np.square(nes * np.cos(thetap))
    alphas = ((np.square(nos) - np.square(nes)) * (np.sin(thetap)) * (np.cos(thetap))) / dens
    betas = (nos * nes) / dens
    gammas = nos / np.sqrt(dens)
    etas = nes * gammas
    nobar = noo(2 * wavep)
    ###etapp = (no*ne)/math.sqrt(den)

    #########For wavefunction(e --> oo) ###################
    funcmat=np.zeros([grdpnt,grdpnt])
    funcmat1=np.zeros([grdpnt,grdpnt])
    for k in np.arange(-5,5):
      for m in np.arange(-5,5):
    #////////////////////////////////////////////////////////
          noi = noo(2*wavep)
          nei = neo(2*wavep)
          thetai= thetap + np.arctan2(-ys + m*dx,-xs + k*dx)
          deni = np.square(noi * np.sin(thetap)) + np.square(nei * np.cos(thetap))
          alphai = ((np.square(noi) - np.square(nei)) * (np.sin(thetap)) * (np.cos(thetap))) / deni
          betai = (noi * nei) / deni
          gammai = noi / np.sqrt(deni)
          etai = nei * gammai
          nobar = noo(2 * wavep)#88888888888888__________________----************************//////////////////////
          funcmat=funcmat+ (np.square(np.absolute(two_photon_wavefunction(xs,ys,-xs+ k*dx,-ys+ m*dx, wp))))*dx*dx
    plt.imshow(funcmat)
    plt.colorbar()
    plt.clim(0.0000001,0.000045) # (0.0000001,0.000045)
    plt.title('Angle =%g' %(angle))
    ##__________________________________________________________________________________________________##
    if not os.path.isdir('static'):
        os.mkdir('static')
    else:
        # Remove old plot files
        for filename in glob.glob(os.path.join('static', '*.png')):
            os.remove(filename)
        # Use time since Jan 1, 1970 in filename in order make
        # a unique filename that the browser has not chached
    plotfile = os.path.join('static', str(time.time()) + '.png')
    plt.savefig(plotfile)
    plt.cla()
    plt.clf()
    return plotfile
#__________________________________________________________________________________________________#
if __name__ == '__main__':
    print SPDC(405, 28.649, 2, 100,388)    #(wavep, angle, L, distz, wp)
