# -*- coding: utf-8 -*-
"""
Programa para verificação de seções de concreto armado submetidas à flexão
simples segundo a NBR-6118:2014, voltado inicialmente à seções retangulares.

@author: GabrielMachado
Last-Modified: 10/07/2022
"""

import numpy as np
import pandas as pd
# import matplotlib.pyplot as plt
# Para resolver o sistema de eq. não lineares
from scipy.optimize import fsolve
#import copy

# Entrada de dados
#

with open('dados.xlsx', 'rb') as target:
    sheet =  pd.read_excel(target, sheet_name='Planilha1')
    data  =  sheet.values
#
nc = np.int_( data[0,0] )         # coluna A: nc
nc1 = nc+1
xc = np.zeros(nc1)
yc = np.zeros(nc1)
xc[0:nc] = data[0:nc,1]           # coluna B: xc
xc[-1] = data[0,1]
yc[0:nc] = data[0:nc,2]           # coluna C: yc
yc[-1] = data[0,2]
As = np.single(data[0,3])         # coluna D: As
ns = np.int_(data[0,4])           # coluna E: ns
xs = np.zeros(ns)
ys = np.zeros(ns)
xs = data[0:ns,5]                 # coluna F: xs
ys = data[0:ns,6]                 # coluna G: ys
rj = data[0:ns,7]                 # coluna H: rj
fck = np.single(data[0,8])        # coluna I: fck
fyk = np.single(data[0,9])        # coluna J: fyk
Es = np.single(data[0,10])        # coluna K: Es
gamac = np.single(data[0,11])     # coluna L: gamac
gamas = np.single(data[0,12])     # coluna M: gamas
Nad = np.single(data[0,13])       # coluna N: Nad
Maxd = np.single(data[0,14])      # coluna O: Maxd
#
#
# Final da entrada de dados

"""
################################################
# dados iniciais de teste
#   geometria da secao
b   = 20
h   = 40
nc = 4
nc1 = nc+1
xc=np.zeros(nc1)
xc[0] = 0
xc[1] = b
xc[2] = b
xc[3] = 0
xc[4] = 0

yc=np.zeros(nc1)
yc[0] = 0
yc[1] = 0
yc[2] = h
yc[3] = h
yc[4] = 0

#   armaduras
rj      = np.array([0.5, 0.5])
As = 2*(np.pi*1.6**2/4)
xs      = np.zeros(2)
xs[0]   = 3.3
xs[1]   = b-3.3

ys=np.zeros(2)
ys[0] = 3.3
ys[1] = 3.3

#   propriedades dos materiais
fck  = 2     #kN/cm²
nbar = 2
Es   = 21000 #kN/cm²
fyk  = 50    #kN/cm²

#   coef. parciais NBR-6118
gamac = 1.4
gamas = 1.15

#   solicitações
Nad  = 0
Maxd = -1.E4

# teste seção t
#del yc, xc
#yc = np.array([0, 0, 10, 15, 95, 110, 120, 120, 110, 95, 15, 10, 0])
#xc = np.array([22.5, 47.5, 47.5, 42.5, 42.5, 70, 70, 0, 0, 27.5, 27.5,
#               22.5, 22.5])
#"""
################################################
# cálculo das const. a1 e a2 para diag. par-retangulo
# verificar se tal procedimento pode ser feito
fcd = fck/gamac
cc  = 210000/Es
fck = fck*cc

if fck<=50:
    epsc2 = 0.002       #def. especificas da nbr 6118
    epscu = 0.0035
    a1    = 1000
    a2    = 250000
else:
    epsc2 = 0.002 + 0.000085*(fck-50)**0.53
    epscu = 0.0026 + 0.035*((90-fck)/100)**4
    nn    = 1.4 + 23.4*((90-fck)/100)**4
    #
    ddx   = epsc2/1000
    x   = 0
    x2  = 0
    x3  = 0
    x4  = 0
    xy  = 0
    x2y = 0
    for i in range(1000):
        y   = 1-(1-x/epsc2)**nn
        x2  = x2 + x*x
        x3  = x3 + x*x*x
        x4  = x4 + x*x*x*x
        xy  = xy + x*y
        x2y = x2y + x*x*y
        x   = x + ddx
    a1 = (x2y-x4*xy/x3)/(x3-x2*x4/x3)
    a2 = -(x2y-x3*a1)/x4
#
xmax = np.amax(xc, axis=0)
xmin = np.amin(xc, axis=0)
ymax = np.amax(yc, axis=0)
ymin = np.amin(yc, axis=0)

# cálculo das propriedades geométricas da seção
Lx = xmax-xmin   #b
Ly = ymax-ymin   #h
Ac = 0
Sx = 0
Jx = 0
Sy = 0
Jy = 0
Jxy = 0
for i in range(nc):
    dx = xc[i+1]-xc[i]
    dy = yc[i+1]-yc[i]
    Ac = Ac + (xc[i]+dx/2)*dy
    Sx = Sx + (xc[i]*(yc[i]+dy/2) + dx*(yc[i]/2+dy/3))*dy
    Jx = Jx + (xc[i]*(yc[i]*(dy+yc[i])+dy*dy/3) + \
               dx*(yc[i]*(yc[i]/2+dy/1.5) + dy*dy/4))*dy
    Sy = Sy + (xc[i]*(xc[i]+dx)+dx*dx/3)*dy/2
    Jy = Jy + (dx**3/4+xc[i]*(dx*dx+xc[i]*(1.5*dx+xc[i])))*dy/3
    Jxy = Jxy + (xc[i]*(xc[i]*(yc[i]+dy/2) + dx*(yc[i]+dy/1.5)) + \
        dx*dx*(yc[i]/3+dy/4))*dy/2
xg   = Sy/Ac
yg   = Sx/Ac
Sxg  = Sx - yg*Ac
Syg  = Sy - xg*Ac
Jxg  = Jx - Ac*yg**2
Jyg  = Jy - Ac*xg**2
Jxyg = Jxy - Ac*xg*yg

# coordenadas no sistema local xg,yg
xc = xc - xg
yc = yc - yg
xs = xs - xg
ys = ys - yg

###################################################

#############################
"""
Função esfor() retorna os esforços resistentes MRxd e NRd, calculados
a partir das funções poli(), reg1(), reg2() - integração da área de 
concreto comprimida -  e aco() - diagrama tensão-deformação das armaduras.
Variáveis:
    X - altura da linha neutra
    Ly - altura da seção
    xc, yc - vértices da seção em relação ao cg.
    ys - coordenadas y das armaduras em relação ao cg.
    epscu, epsc2 - def. específicas do diagrama parábola-retângulo 
            (NBR 6118:2014)
    xg, yg - coord. do centroide
    YS, YI - alturas acima e abaixo da LN
    d - altura útil da seção
    epsS, epsI - def. específicas nos pontos extremos da seção
    b, c - variáveis de ajuste usadas em D0, D1 e D2
    a1, a2 - coeficientes do ajuste por regressão do diagr. parábola/ret
    fcd - resistência à compressão do concreto em projeto
    fyk - resistência à tração característica do aço
    gamas - coeficiente de segurança aço
    rj() - porcentagem de armaduras das barras
    Es - módulo de def. longitudinal do aço
"""
def esfor(As, a1, a2, Es, epscu, epsc2, fcd, fyk, gamas, Ly, rj, X, xc, xg, yc, yg, ys):
    YS = np.abs(np.max(yc))
    YI = np.abs(np.min(ys))
    d  = YS-YI
    # calculo de epsilon S e epsilon I
    if X>(-1E50) and X<=((epscu*d)/(10/1000+epscu)):
        epsS = -10/1000*X/(d-X)
        epsI = 10/1000
    elif X>=((epscu*d)/(10/1000+epscu)) and X<=(d):
        epsS = -epscu
        epsI = epscu*(d-X)/X
    elif X>=(d) and X<=(Ly):
        epsS = -epscu
        epsI = 0
    else:
        epsS = -epsc2*(X/(X-Ly*((epscu-epsc2)/epscu)))
        epsI = -epsc2*((X-Ly)/(X-Ly*((epscu-epsc2)/epscu)))

    b = (epsS-epsI)/(YS-YI)
    c = epsS - b*YS

    xx1, yy1, xx2, yy2 = poli(xc, yc, xg, yg, c, b, epsc2)
    Mrx1, Nr1 = reg1(a1, a2, c, b, fcd, xx1, yy1)
    Mrx2, Nr2 = reg2(fcd, xx2, yy2)
    
    fyd = fyk/gamas
    Mrxas, Nras = aco(As, b, c, Es, fyd, rj, ys)
    
    Mrxd = Mrx1 + Mrx2 + Mrxas
    Nrd = Nr1 + Nr2 + Nras
    
    #print(xx1, yy1)
    #print(Mrx1, Nr1)
    #print(Mrx2, Nr2)
    #print(Mrxas, Nras)
    #print (Mrxd, Nrd)
    #print(epsS, epsI)
    return (Mrxd, Nrd)
#############################

#############################
"""
Função poli() retorna as variáveis xx1,yy1,xx2,yy2 que determinam
as poligonais das regiões 1 e 2 da área comprimida da seção.
Variáveis:
    xc, yc - vértices da seção em relação ao cg.
    xg, yg - coord. do centroide
    c, b - variáveis de ajuste usadas em D0, D1 e D2
    epsc2 - def. específica c2 do diagrama parábola-retângulo 
            (NBR 6118:2014)

Last-Modified: 02/07/2022
Status: Ok - Testado para seções retangulares e T.
"""
def poli(xc, yc, xg, yg, c, b, epsc2):
    y01 = (-c/b)-yg
    y12 = (-epsc2-c)/b-yg

    xx1 = np.copy(xc)
    yy1 = np.copy(yc)

    for i in range(np.size(yy1)):
        if yy1[i]<y01:
            yy1[i] = np.nan
            xx1[i] = np.nan
        elif np.isnan(yy1[i-1])==True:
            yy1[i-1] = y01
            xx1[i-1] = (y01-yc[i-1])/(yc[i]-yc[i-1])*(xc[i]-xc[i-1])+xc[i-1]
        elif i<(np.size(yy1)-1):
            if (yy1[i+1]<y01) and (yy1[i]>y01):
                yy1[i+1] = y01
                xx1[i+1] = (y01-yc[i])/(yc[i+1]-yc[i])*(xc[i+1]-xc[i])+xc[i]
        else:
            continue
    
    for i in range(np.size(yy1)):
        if i<(np.size(yy1)-1):
            if yy1[i]<y12 and yy1[i+1]>y12:
                yy1[i+1] = y12
                xx1[i+1] = (y12-yc[i])/(yc[i+1]-yc[i])*(xc[i+1]-xc[i])+xc[i]
            elif yy1[i]>y12 and yy1[i+1]<y12:
                yy1[i] = y12
                xx1[i] = (y12-yc[i+1])/(yc[i]-yc[i+1])*(xc[i]-xc[i+1])+xc[i+1]
            elif yy1[i]>y12:
                yy1[i] = np.nan
                xx1[i] = np.nan
            else:
                continue
    yy1=yy1[np.logical_not(np.isnan(yy1))]
    xx1=xx1[np.logical_not(np.isnan(xx1))]

    xx2 = np.copy(xc)
    yy2 = np.copy(yc)

    for i in range(np.size(yy2)):
        if yy2[i]<y12:
            yy2[i] = np.nan
            xx2[i] = np.nan
        elif np.isnan(yy2[i-1])==True:
            yy2[i-1] = y12
            xx2[i-1] = (y12-yc[i-1])/(yc[i]-yc[i-1])*(xc[i]-xc[i-1])+xc[i-1]
        elif i<(np.size(yy2)-1):
            if (yy2[i+1]<y12) and (yy2[i]>y12):
                yy2[i+1] = y12
                xx2[i+1] = (y12-yc[i])/(yc[i+1]-yc[i])*(xc[i+1]-xc[i])+xc[i]
        else:
            continue
    yy2=yy2[np.logical_not(np.isnan(yy2))]
    xx2=xx2[np.logical_not(np.isnan(xx2))]
    
    return (xx1, yy1, xx2, yy2)
#############################

#############################
"""
    D0, D1, D2 - coef. geométricos da seção p/ calc. das tensões no concreto
"""
def reg1(a1, a2, c, b, fcd, xx1, yy1):
    D0   = a1*c + a2*c**2
    D1   = a1*b + 2*a2*b*c
    D2   = a2*b**2
    scd  = 0.85*fcd
    Mrx1 = 0
    Nr1  = 0
    for i in range(np.size(xx1)-1):
        x1 = xx1[i]
        y1 = yy1[i]
        x2 = xx1[i+1]
        y2 = yy1[i+1]
        dx = x2 - x1
        dy = y2 - y1
        dy1  = dy/2
        if dy==0:
            continue
        else:
            dy2  = dy*dy
            dy3  = dy2*dy
            G00  = (x1 + dx/2)*dy
            G01  = (x1*(y1+dy1)+dx*(y1/2+dy/3))*dy
            G02  = (x1*(y1*(dy+y1)+dy2/3) + dx*(y1*(y1/2+dy/1.5)+dy2/4))*dy
            G03  = (x1*(y1*(dy2+y1*(1.5*dy+y1))+dy3/4) + 
                    dx*(y1*(0.75*dy2+y1*(dy+y1/2))+dy3/5))*dy
            Mrx1 = Mrx1 + scd*(D0*G01 + D1*G02 + D2*G03)
            Nr1  = Nr1 + scd*(D0*G00 + D1*G01 + D2*G02)
    return (Mrx1, Nr1)
#############################
#############################
def reg2(fcd, xx2, yy2):
    scd  = 0.85*fcd
    Mrx2 = 0
    Nr2  = 0
    for i in range(np.size(xx2)-1):
        x1 = xx2[i]
        y1 = yy2[i]
        x2 = xx2[i+1]
        y2 = yy2[i+1]
        dx = x2 - x1
        dy = y2 - y1
        dy1  = dy/2 
        if dy==0:
            continue
        else:
            G00  = (x1 + dx/2)*dy
            G01  = (x1*(y1+dy1)+dx*(y1/2+dy/3))*dy
            Mrx2 = Mrx2 - scd*G01
            Nr2  = Nr2 - scd*G00
    return (Mrx2, Nr2)
#############################
#############################
def aco(As, b, c, Es, fyd, rj, ys):
    Nras  = 0
    Mrxas = 0
    for i in range(np.size(ys)):
        epsb = b*ys[i] + c
        epsyd = fyd/Es
        if (np.abs(epsb)<=epsyd):
            sig = Es*epsb
        else:
            sig = np.sign(epsb)*fyd
        Nrasi = rj[i]*As*sig
        Nras  = Nras + Nrasi
        Mrxas = Mrxas + Nrasi*ys[i]
    return (Mrxas, Nras)
#############################

def nlsistema(var, *var_aux):
    (X, lamb) = var
    const, epscu, epsc2, rj, xc, xg, yc, yg, ys = var_aux
    As    = const[0]
    a1    = const[1]
    a2    = const[2]
    Es    = const[3]
    fcd   = const[4]
    fyk   = const[5]
    gamas = const[6]
    Ly    = const[7]
    Maxd  = const[8]
    Nad   = const[9]
    Mrxd, Nrd = esfor(As, a1, a2, Es, epscu, epsc2, fcd, fyk, gamas, Ly, rj, X, xc, xg, yc, yg, ys)
    f = lamb*(Mrxd) - Maxd
    g = lamb*(Nrd) - Nad
    return[f,g]
const = np.array([As, a1, a2, Es, fcd, fyk, gamas, Ly, Maxd, Nad])
lambi = 0.5 #lambda inicial
Xi = 0 #altura da LN inicial
s0 = np.array([Xi, lambi])
var_aux = (const, epscu, epsc2, rj, xc, xg, yc, yg, ys)
s  = fsolve(nlsistema, s0, var_aux)
print(s)
