#%%

# INFOS GÉNÉRALES
# Choisir dans les choix le calcul désiré (contraintes normales ou de cisaillement)
# Les calculs pour DCL et forces internes se font peut importe le choix

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

# Données générales

#Choix du calcul 
# 1 = Contraintes normales
# 2 = Contraintes de cisaillement
choix = 1

step = 0.001
x = np.arange(0,1.51+step,step)

#Fonctions à intégrer
def wx(x):
    return (-266.01*(x**6)) + (1099.1*(x**5)) - (1707.7*(x**4)) + (1227.1*(x**3)) - (413.43*(x**2)) + (41.818*x) + 63.532

def x_wx(x):
    return wx(x)*x

#------------------------DCL --------------------------

#Calcul de la force equivalente pour le DCL
centroide = (1/(integrate.quad(wx,0,1.51)[0]))*integrate.quad(x_wx,0,1.51)[0]
force_equivalente = integrate.quad(wx,0,1.51)[0]

# Somme des forces
Ay = force_equivalente
M0 = force_equivalente*centroide

#-----------------------Forces internes-------------------

# V(x)
Vx = []

for i in x:
    Vx.append(-74.0678 + integrate.quad(wx,0,i)[0])

figure1, figure_v = plt.subplots()
figure_v.plot(x,Vx)
figure_v.set(title = "Effort tranchant selon la distance sur l'aile", xlabel="Distance x sur l'aile (m)", ylabel = 'Valeur de V(x) (N)')

# M(x)
Mx = []

for i in x:
    if i != 0:
        centroide_2 = (1/(integrate.quad(wx,0,i)[0]))*integrate.quad(x_wx,0,i)[0]
        value = (-Ay*i) + M0 + (i-centroide_2)*(integrate.quad(wx,0,i)[0])
        Mx.append (value)
    else:
        Mx.append(47.6182)

figure2, figure_m = plt.subplots()
figure_m.plot(x,Mx)
figure_m.set(title = "Moment selon la distance sur l'aile", xlabel="Distance x sur l'aile (m)", ylabel = 'Valeur de M(x) (N*m)')

#--------------------------Contraintes normales (compression et traction)---------
if ( choix ==1):
    # Contraintes normales 

    # Variables necessaires
    t = 3.5*10**-4

    def c(x): #Longueur de la corde au point x
        return -0.1762*x+0.606

    def Iz(x): #Valeur de Ix dans les donnees
        return 0.0051752 *(c(x)**3)*t

    def yc(x): #Centroide en y, utilise dans Excel seulement
        return 0.046128*c(x)

    #Valeurs trouvees avec Excel
    def distance_max_sup(x):
        return 0.08138347 * c(x) # Valeur a x/c = 0.38885365

    def distance_max_inf(x): # Valeur a x/c = 0.07820223
        return -0.06393017 * c(x)

    # Compression
    contrainte_max_compression = []

    loop_nb = 0
    for i in x:
        value = (-(Mx[loop_nb])*distance_max_sup(i)) / (Iz(i)) 
        contrainte_max_compression.append(value)
        loop_nb = loop_nb + 1

    figure3, figure_c1 = plt.subplots()
    figure_c1.plot(x,contrainte_max_compression)
    figure_c1.set(title = "Contrainte en compression max selon la distance sur l'aile", xlabel="Distance x sur l'aile (m)", ylabel = 'Valeur de la contrainte en compression (Pa)')

    #Compression

    contrainte_max_traction = []

    loop_nb = 0
    for i in x:
        value = (-(Mx[loop_nb])*distance_max_inf(i)) / (Iz(i)) 
        contrainte_max_traction.append(value)
        loop_nb = loop_nb + 1

    figure4, figure_c2 = plt.subplots()
    figure_c2.plot(x,contrainte_max_traction)
    figure_c2.set(title = "Contrainte en traction max selon la distance sur l'aile", xlabel="Distance x sur l'aile (m)", ylabel = 'Valeur de la contrainte en traction (Pa)')

    #Valeurs finales

    Contrainte_compression = min(contrainte_max_compression)
    Contrainte_traction = max(contrainte_max_traction)

    print ("Contrainte en compression maximale: ", format(Contrainte_compression, "e"), " Pa" )
    print ("Contrainte en traction maximale: ", format(Contrainte_traction, "e"), " Pa")
    print("Les 2 valeurs sont maximales au point x=0")

    #Graphique des contraintes selon la distance y sur le profil de l'aile
    #Pour prouver que la variation des contraintes est linéaire

    #Superieur
    step = 0.001
    y_sup = np.arange(0,distance_max_sup(0)+step, step)

    contrainte_compression_selon_y = []

    loop_nb = 0
    for i in range(0,len(y_sup)):
        value = (-(Mx[i])*y_sup[i]) / (Iz(0)) 
        contrainte_compression_selon_y.append(-value)
        #Le signe négatif a été rajouté juste pour que le graphique soit plus concis
        loop_nb = loop_nb + 1

    figure5, figure_c3 = plt.subplots()
    figure_c3.plot(contrainte_compression_selon_y, y_sup, label="compression", color="red")
    figure_c3.set(ylabel="Distance y sur le profil (m)", xlabel = 'Contrainte (Pa)')

    #Inferieur
    step = 0.001
    y_inf = np.arange(distance_max_inf(0), 0+step, step)

    contrainte_traction_selon_y = []

    loop_nb = 0
    for i in range(0,len(y_inf)):
        value = (-(Mx[i])*y_inf[i]) / (Iz(0)) 
        contrainte_traction_selon_y.append(value)
        loop_nb = loop_nb + 1

    figure_c3.plot(contrainte_traction_selon_y, y_inf, label ="traction", color = "blue")
    figure_c3.legend(loc="upper left")


    plt.show()

#---------------------------------------Contraintes de cisaillement-------------
elif (choix == 2):
#Contraintes de cisaillement

    import csv

    #Import CSV profile data
    datafile = open('profil.csv')
    myreader = csv.reader(datafile, delimiter=';')
    next(myreader, None)
    lst = list(myreader)
    arr_txt=np.array(lst)
    arr = np.asarray(arr_txt,dtype=float)

    #Separate x and y values
    profile_x = arr[:,0]
    profile_y = arr[:,1]

    #Show profile
    figure3, figure_profile = plt.subplots()
    figure_profile.scatter(profile_x,profile_y)
    figure_profile.set(title = "Profil de l'aile", xlabel="x/c", ylabel = 'y/c')    

    # Variables necessaires
    #t = 3.5*10**-4

    def c(x): #Longueur de la corde au point x
        return -0.1762*x+0.606

    def Iz(x): #Valeur de Ix dans les donnees
        return 9.8490*(10**-5) *(c(x)**4)
    
    def yc(x): #Centroide en y
        return 0.053017*c(x)
    
    # ---------------------------------------Calcul de Q-------------------------------
    step_2 = 0.001
    value_y =  0.053017 #centroide general profil

    profile_top = []

    #Points pour le profil de la section plane
    for i in range(len(arr)):
        if (arr[i,1] >= value_y):
            profile_top.append([arr[i,0],arr[i,1]])

    profile_top_arr = np.array(profile_top)

    profile_2_x = profile_top_arr[:,0]
    profile_2_y = profile_top_arr[:,1]

    #Formule pour le profil au dessus du centroide
    #La formule est cherchee a 90 degre (donc x(y)) pour trouver le centroide x (qui est y au final)
    def formula_profile_top(x):
        N=10 #Degre pour la formule
        value = 0
        coeffs = np.polyfit(profile_2_y, profile_2_x, N)
        for i in range(N+1):
            value = value + coeffs[i]*(x**N)
            N = N-1
        value = value - value_y # L'aire est entre la courbe trouvee, et une droite passant par value_y (le centroide sans angle)
        return (value)
    
    def formula_profile_top_y(x):
        return (formula_profile_top(x)*x)
    
    #Calcul du centroide en y 
    Aire = ((integrate.quad(formula_profile_top,min(profile_2_y), max(profile_2_y))[0]))
    centroide_y = ( (1/Aire)*integrate.quad(formula_profile_top_y,min(profile_2_y), max(profile_2_y))[0])

    # Affichage du profil au dessus du centroide
    figure4, figure_pr = plt.subplots()
    figure_pr.scatter(profile_2_x,profile_2_y)
    figure_pr.set(title = "Section coupée du profil de l'aile", xlabel="x/c", ylabel = 'y/c')  

    #Position centroide par rapport a l'axe neutre (en 1/m)
    position = centroide_y - value_y

    #Calcul final Q (en 1/m^3. Conversion a la fin avec valeur de c(x))
    Q = position*Aire

    #-------------------------------Calcul de t (largeur de la section) ---------------------------------------------------------
    #Ici, t est egal a la distance y a l'axe neutre
    # Largeur_t est en x/c, Conversion a la fin avec c(x)
    largeur_t = max(profile_2_x) - min(profile_2_x)

    #Calcul contrainte cisaillement
    contrainte_max_cisaillement = []
    loop_nb = 0
    for i in x:
        value = (Vx[loop_nb]*Q*(c(i)**3))/(Iz(i)*largeur_t*c(i))
        contrainte_max_cisaillement.append(value)
        loop_nb = loop_nb + 1

    #Affichage graphique contrainte cisaillement selon x
    figure5, figure_cis = plt.subplots()
    figure_cis.scatter(x,contrainte_max_cisaillement)
    figure_cis.set(title = "Contrainte de cisaillement selon la position sur l'aile sur l'axe neutre", xlabel="distance x sur l'aile (m)", ylabel = 'contrainte cisaillement max (Pa)')  

    plt.show()

    #Valeur de la contrainte maximale finale
    print ("Contraine de cisaillement max: ", min(contrainte_max_cisaillement)/1000, " KPa")


# %%
