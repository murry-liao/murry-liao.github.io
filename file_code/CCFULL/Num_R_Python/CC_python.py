# -*- encoding: utf-8 -*-

'''
@File    :   CC_python.py
@Time    :   2025/01/08 10:20:35
@Author  :   Liao ZeHong, Kouichi Hagino
@Version :   1.0
@Contact :   liaozh26@mail2.sysu.edu.cn
'''
# here put the import lib
import time
import numpy as np
from mpmath import coulombf, coulombg
from scipy.signal import argrelmax, argrelmin
from sympy.physics.quantum.cg import CG
import matplotlib.pyplot as plt
from math import pi, sqrt
import pandas as pd

def WS_parameter_Aw():
    gama = 0.95*(1 - 1.8 * (Apro-2*Zpro)/Apro* (Atar-2*Ztar)/Atar)
    R1 = 1.2 * Apro**(1./3) - 0.09
    R2 = 1.2 * Atar**(1./3) - 0.09
    a = 1. / (1.17 * (1 + 0.53 * (Apro**(-1./3) + Atar**(-1./3))))
    R_hbar = R1 * R2 / (R1 + R2)
    V0 = 16 * 3.1415926 * gama * R_hbar * a
    R0 = R1 + R2

    return V0, a, R0

def WS_parameter():
    return V0, A0, R0

def Vn(r):
    V0, A0, R0 = WS_parameter()
    R12 = R0 * (Apro**(1./3)+Atar**(1./3))
    return  -V0 / (1 + np.exp((r - R12) / A0))
    
def dVn_dr(r):
    V0, A0, R0 = WS_parameter()
    R12 = R0 * (Apro**(1./3)+Atar**(1./3))
    return V0/A0 * np.exp((r - R12) / A0) / (1 + np.exp((r - R12) / A0))**2

def VnCC(r, Xt):
    V0, A0, R0 = WS_parameter()
    R12 = R0 * (Apro**(1./3)+Atar**(1./3))
    return -V0 / (1 + np.exp((r - R12 - Xt) / A0))

def Vc(r):
    return Zpro * Ztar / r * Hbar / 137.0        #Zpro * Ztar * e2 / r

def Vr(r, l):
    return h2m * (l * (l + 1.)) / r ** 2

def V(r, l):
    return Vn(r) + Vc(r) + Vr(r, l)

def grotation():
    global erott, erotp

    for it in range(Ntar + 1):
        for jt in range(Ntar + 1):
            if it == jt - 1 or it == jt + 1 or it == jt:
                bett[it][jt] = Beta2T
        erott[it] = 2.0 * it * (2.0 * it + 1.0) / 6.0 * E2T
        if IVIBROTT == 1:
          print(f"target rotation exciation of {2*it} + = :" + str(erott[it]) + " MeV")

    if IVIBROTT == 1 and Ntar > 1:
        ans = input("Generalised E2 couplings for the target rotor (y/n)? ").strip().lower()
        if ans == 'y':
            print("\n**** Generalised couplings for the target rotor ****")
            print(f"Beta2Targ = {Beta2T}")

            for it in range(1, Ntar):
                jt = it + 1
                ans = input(f"Modify beta2 for the transition from {2 * it}+ to {2 * jt}+ (y/n)? ").strip().lower()
                if ans == 'y':
                    bett[it][jt] = float(input("   BETA2=? "))
                    bett[jt][it] = bett[it][jt]
                    print(f"Beta2 for the transition from {2 * it}+ to {2 * jt}+ = {bett[it][jt]}")

            print("\nReorientation terms:")
            for it in range(0, Ntar):
                jt = it + 1
                ans = input(f"Modify beta2 for the transition from {2 * it}+ to {2 * jt}+ (y/n)? ").strip().lower()
                if ans == 'y':
                    bett[it][jt] = float(input("   BETA2=? "))
                    bett[jt][it] = bett[it][jt]
                    print(f"Beta2 for the transition from {2 * it}+ to {2 * jt}+ = {bett[it][jt]}")

            print("\nExcitation energies:")
            for i in range(2, Ntar + 1):
                print(f"    Energy of the {2 * i}+ state for the pure rotor = {erott[i]}")
                ans = input("   Modify this energy (y/n)? ").strip().lower()
                if ans == 'y':
                    erott0 = erott[i]
                    erott[i] = float(input("   Energy=? "))
                    print(f"Energy of the {2 * i}+ state = {erott[i]} (rotor: {erott0})")

    for ip in range(Npro + 1):
        for jp in range(Npro + 1):
            if ip == jp - 1 or ip == jp + 1 or ip == jp:
                betp[ip][jp] = Beta2P
        erotp[ip] = 2.0 * ip * (2.0 * ip + 1.0) / 6.0 * E2P
        if IVIBROTP == 1:
          print(f"projectile rotation exciation of {2*ip} + = :" + str(erotp[ip]) + " MeV")

    if IVIBROTP == 1 and Npro > 1:
        ans = input("Generalised E2 couplings for the projectile rotor (y/n)? ").strip().lower()
        if ans == 'y':
            print("\n**** Generalised couplings for the projectile rotor ****")
            print(f"Beta2Proj = {Beta2P}")

            for ip in range(1, Npro):
                jp = ip + 1
                ans = input(f"Modify beta2 for the transition from {2 * ip}+ to {2 * jp}+ (y/n)? ").strip().lower()
                if ans == 'y':
                    betp[ip][jp] = float(input("   BETA2=? "))
                    betp[jp][ip] = betp[ip][jp]
                    print(f"Beta2 for the transition from {2 * ip}+ to {2 * jp}+ = {betp[ip][jp]}")

            print("\nReorientation terms:")
            for ip in range(1, Npro + 1):
                jp = ip
                ans = input(f"Modify beta2 for the transition from {2 * ip}+ to {2 * jp}+ (y/n)? ").strip().lower()
                if ans == 'y':
                    betp[ip][jp] = float(input("   BETA2=? "))
                    betp[jp][ip] = betp[ip][jp]
                    print(f"Beta2 for the transition from {2 * ip}+ to {2 * jp}+ = {betp[ip][jp]}")

            print("\nExcitation energies:")
            for i in range(2, Npro + 1):
                print(f"    Energy of the {2 * i}+ state for the pure rotor = {erotp[i]}")
                ans = input("   Modify this energy (y/n)? ").strip().lower()
                if ans == 'y':
                    erotp0 = erotp[i]
                    erotp[i] = float(input("   Energy=? "))
                    print(f"Energy of the {2 * i}+ state = {erotp[i]} (rotor: {erotp0})")

    return

def anharmonicity():
    global betcahv, betcahv, omeahv
    global betcahv2, betcahv2, omeahv2
    global betcahvp, betcahv, omeahvp
    """

    Returns:
      betnahv   (list):  Anharmonic coupling parameters for the first mode of the target nucleus.
      betcahv   (list):  Anharmonic coupling parameters for the first mode of the target nucleus.
      omeahv    (list):  Excitation energy for the first mode of the target nucleus.
      betnahv2  (list):  Anharmonic coupling parameters for the second mode of the target nucleus.
      betcahv2  (list):  Anharmonic coupling parameters for the second mode of the target nucleus.
      omeahv2   (list):  Excitation energy for the second mode of the target nucleus.
      betnahvp  (list):  Anharmonic coupling parameters for the projectile nucleus.
      betcahvp  (list):  Anharmonic coupling parameters for the projectile nucleus.
      omeahvp   (list):  Excitation energy for the projectile nucleus.
    """

    for it in range(Ntar + 1):
        for jt in range(Ntar + 1):
            if it == jt - 1:
                betnahv[it][jt] = BetaTn * (jt ** 0.5)
                betcahv[it][jt] = BetaT * (jt ** 0.5)
            elif jt == it - 1:
                betnahv[it][jt] = BetaTn * (it ** 0.5)
                betcahv[it][jt] = BetaT * (it ** 0.5)
        omeahv[it] = it * OmegaT

    if IVIBROTT == 0 and Ntar > 1:
        ans = input("AHV couplings for the first mode in the target phonon (y/n)? ").strip().lower()
        if ans == 'y':
            print("\n**** AHV Couplings in the target (the 1st mode) ****")
            sign = '+' if (-1) ** LambdaT == 1 else '-'
            for it in range(Ntar + 1):
                for jt in range(Ntar + 1):
                    if it > jt or (it == 0 and jt == 0) or (it == 0 and jt == 1):
                        continue
                    print(f"\nTransition from the {LambdaT}{sign}^{it} to the {LambdaT}{sign}^{jt} state:")
                    print(f"   beta_N and beta_C in the HO limit: {betnahv[it][jt]}, {betcahv[it][jt]}")
                    ans = input("    Modify these beta_N and/or beta_C (y/n)? ").strip().lower()
                    if ans == 'y':
                        betnahv[it][jt], betcahv[it][jt] = map(float, input("   beta_N and beta_C =? ").split())


            print("\nExcitation energy for the first mode:")
            for i in range(2, Ntar + 1):
                print(f"    Energy of the {i}-phonon state in the HO: {i * OmegaT}")
                ans = input("   Modify this energy (y/n)? ").strip().lower()
                if ans == 'y':
                    omeahv[i] = float(input("   Energy=? "))

    for it in range(NphononT2 + 1):
        for jt in range(NphononT2 + 1):
            if it == jt - 1:
                betnahv2[it][jt] = BetaT2n * (jt ** 0.5)
                betcahv2[it][jt] = BetaT2 * (jt ** 0.5)
            elif jt == it - 1:
                betnahv2[it][jt] = BetaT2n * (it ** 0.5)
                betcahv2[it][jt] = BetaT2 * (it ** 0.5)
        omeahv2[it] = it * OmegaT2

    if NphononT2 > 1:
        ans = input("AHV couplings for the second mode in the target phonon (y/n)? ").strip().lower()
        if ans == 'y':
            print("\n**** AHV Couplings in the target (the 2nd mode) ****")
            sign = '+' if (-1) ** LambdaT2 == 1 else '-'
            for it in range(NphononT2 + 1):
                for jt in range(NphononT2 + 1):
                    if it > jt or (it == 0 and jt == 0) or (it == 0 and jt == 1):
                        continue
                    print(f"\nTransition from the {LambdaT2}{sign}^{it} to the {LambdaT2}{sign}^{jt} state:")
                    print(f"   beta_N and beta_C in the HO limit: {betnahv2[it][jt]}, {betcahv2[it][jt]}")
                    ans = input("    Modify these beta_N and/or beta_C (y/n)? ").strip().lower()
                    if ans == 'y':
                        betnahv2[it][jt], betcahv2[it][jt] = map(float, input("   beta_N and beta_C =? ").split())


            print("\nExcitation energy for the second mode:")
            for i in range(2, NphononT2 + 1):
                print(f"    Energy of the {i}-phonon state in the HO: {i * OmegaT2}")
                ans = input("   Modify this energy (y/n)? ").strip().lower()
                if ans == 'y':
                    omeahv2[i] = float(input("   Energy=? "))


    for it in range(Npro + 1):
        for jt in range(Npro + 1):
            if it == jt - 1:
                betnahvp[it][jt] = BetaPn * (jt ** 0.5)
                betcahvp[it][jt] = BetaP * (jt ** 0.5)
            elif jt == it - 1:
                betnahvp[it][jt] = BetaPn * (it ** 0.5)
                betcahvp[it][jt] = BetaP * (it ** 0.5)
        omeahvp[it] = it * OmegaP

    if IVIBROTP == 0 and Npro > 1:
        ans = input("AHV couplings for the projectile phonon (y/n)? ").strip().lower()
        if ans == 'y':
            print("\n**** AHV Couplings in the projectile ****")
            sign = '+' if (-1) ** LambdaP == 1 else '-'
            for it in range(Npro + 1):
                for jt in range(Npro + 1):
                    if it > jt or (it == 0 and jt == 0) or (it == 0 and jt == 1):
                        continue
                    print(f"\nTransition from the {LambdaP}{sign}^{it} to the {LambdaP}{sign}^{jt} state:")
                    print(f"   beta_N and beta_C in the HO limit: {betnahvp[it][jt]}, {betcahvp[it][jt]}")
                    ans = input("    Modify these beta_N and/or beta_C (y/n)? ").strip().lower()
                    if ans == 'y':
                        betnahvp[it][jt], betcahvp[it][jt] = map(float, input("   beta_N and beta_C =? ").split())


            print("\nExcitation energy for the projectile phonon:")
            for i in range(2, Npro + 1):
                print(f"    Energy of the {i}-phonon state in the HO: {i * OmegaP}")
                ans = input("   Modify this energy (y/n)? ").strip().lower()
                if ans == 'y':
                    omeahvp[i] = float(input("   Energy=? "))

    return 

def mutual():
    global Nlevel, imutual

    print("\nMutual excitations in the *target* nucleus")
    ans = input("Include the mutual excitations (y/n)? ").strip().lower()

    if ans == 'n':

        imut = 0
        Nlevel = (Ntar + NphononT2 + 1) * (Npro + 1)
        print("No mutual excitations in the target are included.")
        for i in range(Ntar + 1):
            for j in range(NphononT2 + 1):
                if i == 0 or j == 0:
                    imutual[i][j] = 1
                    
    else:

        ans = input("All the possible mutual excitation channels (y/n)? ").strip().lower()
        if ans == 'y':

            imut = 1
            Nlevel = (Ntar + 1) * (NphononT2 + 1) * (Npro + 1)
            print("All the possible mutual excitation channels in the target are included.")
            for i in range(Ntar + 1):
                for j in range(NphononT2 + 1):
                    imutual[i][j] = 1
        else:
            # 手动选择互激态
            imut = 2
            sign1 = '+' if (-1) ** LambdaT == 1 else '-'
            sign2 = '+' if (-1) ** LambdaT2 == 1 else '-'
            for i in range(Ntar + 1):
                for j in range(NphononT2 + 1):
                    if i == 0 or j == 0:
                        imutual[i][j] = 1
                        continue
                    ans = input(f"Include ({LambdaT}{sign1}^{i}, {LambdaT2}{sign2}^{j}) state? (y/n) ").strip().lower()
                    if ans == 'n':
                        imutual[i][j] = 0
                    else:
                        imutual[i][j] = 1

    Nlevel = 0
    print("\nExcited states in the target to be included:")
    sign1 = '+' if (-1) ** LambdaT == 1 else '-'
    sign2 = '+' if (-1) ** LambdaT2 == 1 else '-'
    for i in range(Ntar + 1):
        for j in range(NphononT2 + 1):
            if imutual[i][j] == 1:
                print(f"    ({LambdaT}{sign1}^{i}, {LambdaT2}{sign2}^{j}) state")
                Nlevel = Nlevel + 1

    

    Nlevel = Nlevel* (Npro + 1)
    print('Nlevel = ', Nlevel)

    df = pd.DataFrame(imutual)
    print('imutual:', file=open(Filename_output, "a"))
    print(df[0:6][0:6], file=open(Filename_output, "a"))


    return

def coupled_matrix0():
  global Ev, Evec
  # Initialize coupled_matrix
  # It should be noted that in Python, indexing starts from 0, 
  # so A[0][0] represents the first row and first column of the matrix.
  print("preparation for the nuclear coupling matrix")
  A = np.zeros((Nlevelmax, Nlevelmax), dtype=float)

  Ndim = Nlevel - Ntrans
  print(f"Ndim = {Ndim}")
  if Ndim == 1:
      return
  i = -1
  for ip in range(0, Npro + 1):
      for it in range(0, Ntar + 1):
          for it2 in range(0, NphononT2 + 1):
              if  NphononT2 != 0 and imutual[it][it2] == 0:
                  continue
              i = i + 1
              j = -1
              for jp in range(0, Npro + 1):
                  for jt in range(0, Ntar + 1):
                      for jt2 in range(0, NphononT2 + 1):
                          if  NphononT2 != 0 and imutual[jt][jt2] == 0:
                              continue
                          j = j + 1
                          if(i > j):
                              A[i][j] = A[j][i]
                              continue
                          # Calculate coupling term C
                          C = 0.
                          if  ip == jp and it2 == jt2:
                              if IVIBROTT == 0:
                                  C = Rtar * betnahv[it][jt]/np.sqrt(4.0*pi)
                              else:
                                  C = Rtar * bett[it][jt] * sqrt((2 * 2 * it + 1) * 5 * (2 * 2 * jt + 1) /4. / pi)\
                                      * CG(2 * it, 0, 2, 0, 2 * jt, 0).doit()**2 / (2 * 2 * jt + 1)
                                  
                                  C+= Rtar * Beta4T * sqrt((2 * 2 * it + 1) * 9 * (2 * 2 * jt + 1) /4. / pi)\
                                      * CG(2 * it, 0, 4, 0, 2 * jt, 0).doit()**2 / (2 * 2 * jt + 1)
                          if  ip == jp and it == jt:
                              C = C + Rtar * betnahv2[it2][jt2]/np.sqrt(4.0*pi)
                          if  it == jt and it2 == jt2:
                              if IVIBROTP == 0:
                                  C = C + Rpro * betnahvp[ip][jp]/np.sqrt(4.0*pi)
                              else:
                                  C += Rpro * betp[ip][jp] * sqrt((2 * 2 * ip + 1) * 5 * (2 * 2 * jp +1) /4. /pi)\
                                      *CG(2 * ip, 0, 2, 0, 2 * jp, 0).doit()**2 / (2 * 2 * jp + 1)
                                  C += Rpro * Beta4P * sqrt((2 * 2 * ip + 1) * 9 * (2 * 2 * jp +1) /4. /pi)\
                                      *CG(2 * ip, 0, 4, 0, 2 * jp, 0).doit()**2 / (2 * 2 * jp + 1)
                          
                          A[i][j] = C
                          
  # 转换为 DataFrame  打印 DataFrame
  df = pd.DataFrame(A)
  print('coupled_matrix0:', file=open(Filename_output, "a"))
  print(df.iloc[:Ndim, :Ndim], file=open(Filename_output, "a"))

  Ev, Evec = np.linalg.eig(A)
  print('englevalue:', Ev[0:Ndim], file=open(Filename_output, "a"))
  print('englevector:', file=open(Filename_output, "a"))
  print(Evec[0:Ndim][0:Ndim], file=open(Filename_output, "a"))

def Fct(r):
  """_summary_
      Coulomb coupling form factor for the target phonon excitation (first mode)
  Args:
      r (_type_): _description_
  """
  Lambda = LambdaT
  result = 0.
  if r > Rtar:
      result = 3./(2*Lambda + 1) * Zpro * Ztar / 137. *Hbar /r * (Rtar/r)**Lambda
  else:
      result = 3./(2*Lambda + 1) * Zpro * Ztar / 137. *Hbar /Rtar * (r/Rtar)**Lambda
  
  result = result * BetaT / sqrt(4. * pi)

  return result

def Fct2(r):
  """_summary_
      Coulomb coupling form factor for target excitation (rotational E2 coupling)
  Args:
      r (_type_): _description_
  """
  Lambda = 2
  result = 0.
  if r > Rtar:
      result = 3./(2*Lambda + 1) * Zpro * Ztar / 137. *Hbar /r * (Rtar/r)**Lambda
  else:
      result = 3./(2*Lambda + 1) * Zpro * Ztar / 137. *Hbar /Rtar * (r/Rtar)**Lambda
  
  #result = result * (Beta2T + 2. * sqrt(5. / pi) * Beta2T**2 / 7.)

  return result

def Fct3(r):
  """_summary_
    Coulomb coupling form factor for target excitation (rotational E3 coupling)
  Args:
      r (_type_): _description_
  """
  Lambda = 3
  result = 0.
  if r > Rtar:
      result = 3./(2*Lambda + 1) * Zpro * Ztar / 137. *Hbar /r * (Rtar/r)**Lambda
  else:
      result = 3./(2*Lambda + 1) * Zpro * Ztar / 137. *Hbar /Rtar * (r/Rtar)**Lambda
  
  result = result * (Beta4T + 7. * Beta2T**2 /7./sqrt(pi))

  return result

def Fct4(r):
  """_summary_
    Coulomb coupling form factor for target excitation (rotational E4 coupling)
  Args:
      r (_type_): _description_
  """
  Lambda = 4
  result = 0.
  if r > Rtar:
      result = 3./(2*Lambda + 1) * Zpro * Ztar / 137. *Hbar /r * (Rtar/r)**Lambda
  else:
      result = 3./(2*Lambda + 1) * Zpro * Ztar / 137. *Hbar /Rtar * (r/Rtar)**Lambda
  
  result = result * (Beta4T + 9. * Beta2T**2 /7./sqrt(pi))

  return result

def Fcp(r):
  """_summary_
    Coulomb coupling form factor for projectile phonon excitation
  Args:
      r (_type_): _description_
  """
  Lambda = LambdaP
  result = 0.
  if r > Rpro:
      result = 3./(2*Lambda + 1) * Zpro * Ztar / 137. *Hbar /r * (Rpro/r)**Lambda
  else:
      result = 3./(2*Lambda + 1) * Zpro * Ztar / 137. *Hbar /Rpro * (r/Rpro)**Lambda
  
  result = result * BetaP / sqrt(4. * pi)

  return result

def Fcp2(r):
  """_summary_
      Coulomb coupling form factor for projectile excitation (rotational E2 coupling)
  Args:
      r (_type_): _description_
  """
  Lambda = 2
  result = 0.
  if r > Rpro:
      result = 3./(2*Lambda + 1) * Zpro * Ztar / 137. *Hbar /r * (Rpro/r)**Lambda
  else:
      result = 3./(2*Lambda + 1) * Zpro * Ztar / 137. *Hbar /Rpro * (r/Rpro)**Lambda
  
  #result = result * (Beta2P + 2. * sqrt(5. / pi) * Beta2P**2 / 7.)

  return result

def Fcp3(r):
  """_summary_
    Coulomb coupling form factor for projectile excitation (rotational E3 coupling)
  Args:
      r (_type_): _description_
  """
  Lambda = 3
  result = 0.
  if r > Rpro:
      result = 3./(2*Lambda + 1) * Zpro * Ztar / 137. *Hbar /r * (Rpro/r)**Lambda
  else:
      result = 3./(2*Lambda + 1) * Zpro * Ztar / 137. *Hbar /Rpro * (r/Rpro)**Lambda
  
  result = result * (Beta4P + 7. * Beta2P**2 /7./sqrt(pi))

  return result

def Fcp4(r):
  """_summary_
    Coulomb coupling form factor for projectile excitation (rotational E4 coupling)
  Args:
      r (_type_): _description_
  """
  Lambda = 4
  result = 0.
  if r > Rpro:
      result = 3./(2*Lambda + 1) * Zpro * Ztar / 137. *Hbar /r * (Rpro/r)**Lambda
  else:
      result = 3./(2*Lambda + 1) * Zpro * Ztar / 137. *Hbar /Rpro * (r/Rpro)**Lambda
  
  result = result * (Beta4P + 9. * Beta2P**2 /7./sqrt(pi))

  return result

def Fctt(r):
    """_summary_
      Coulomb coupling form factor for the second target phonon excitation
    Args:
        r (_type_): _description_
    """
    Lambda = LambdaT2
    result = 0.
    if r > Rtar:
        result = 3./(2*Lambda + 1) * Zpro * Ztar / 137. *Hbar /r * (Rtar/r)**Lambda
    else:
        result = 3./(2*Lambda + 1) * Zpro * Ztar / 137. *Hbar /Rtar * (r/Rtar)**Lambda
    
    result = result * BetaT2 / sqrt(4. * pi)

    return result

def Fct2v(r):
    """_summary_
      coulomb coupling form factor for target excitation
    Args:
        r (_type_): _description_
    """
    Lambda = 2
    result = 0.
    if r > Rtar:
        result = 3./(2*Lambda + 1) * Zpro * Ztar / 137. * Hbar /r * (Rtar/r)**Lambda
    else:
        result = 3./(2*Lambda + 1) * Zpro * Ztar / 137. * Hbar /Rtar * (r/Rtar)**Lambda
    
    return result

def Fcp2v(r):
    """_summary_
      coulomb coupling form factor for projectile excitation
    Args:
        r (_type_): _description_
    """
    Lambda = 2
    result = 0.
    if r > Rpro:
        result = 3./(2*Lambda + 1) * Zpro * Ztar / 137. * Hbar /r * (Rpro/r)**Lambda
    else:
        result = 3./(2*Lambda + 1) * Zpro * Ztar / 137. * Hbar /Rtar * (r/Rpro)**Lambda
    
    return result

def Ftrans(r):
  """_summary_
  Coupling form factor for transfer reactions
  Args:
      r (_type_): _description_
  """
  result = 0.
  result = Ftr * dVn_dr(r)
  return result

def CoupledMatrix(r):
    """_summary_
        Coupling matrix
    Args:
        r (_type_): _description_
        cpot (_type_): _description_
    """
    
    cpot = np.zeros((Nlevelmax, Nlevelmax), dtype=float)
    Ndim = Nlevel - Ntrans
    #print(f"Npro = {Npro}, Ntar = {Ntar}, NphononT2 = {NphononT2}")
    if Nlevel == 1:
        return
    if Ndim != 1:
        i = -1
        for ip in range(0, Npro + 1):
            for it in range(0, Ntar + 1):
                for it2 in range(0, NphononT2 + 1):
                    if  NphononT2 != 0 and imutual[it][it2] == 0:
                        continue
                    i = i + 1
                    j = -1
                    for jp in range(0, Npro + 1):
                        for jt in range(0, Ntar + 1):
                            for jt2 in range(0, NphononT2 + 1):
                                if  NphononT2 != 0 and imutual[jt][jt2] == 0:
                                    continue
                                j = j + 1
                                if(i > j):
                                    cpot[i][j] = cpot[j][i]
                                    continue
                                C = 0.
                                # Nuclear coupling 
                                for k in range(0, Ndim):
                                    C = C + VnCC(r, Ev[k]) * Evec[i][k] * Evec[j][k]
                                    
                                # Coulomb coupling 
                                # target
                                if  ip == jp and it2 == jt2:
                                    if IVIBROTT == 0:
                                        if it != jt:
                                            cc = betcahv[it][jt] * Fct(r)
                                            if BetaT != 0.:
                                                cc = cc / BetaT
                                        else:
                                            cc = betcahvp[it][jt] * Fct2v(r) / sqrt(4. * pi)

                                        C = C + cc
                                    else:
                                        C = C + sqrt((2 * 2 * it + 1) * 5 * (2 * 2 * jt + 1) /4. / pi)\
                                          * CG(2 * it, 0, 2, 0, 2 * jt, 0).doit()**2 / (2 * 2 * jt + 1) * Fct2(r)\
                                          * (bett[it][jt] + 2. * sqrt(5./pi) * bett[it][jt]**2./7.)
                                        
                                        C = C + sqrt((2 * 2 * it + 1) * 9 * (2 * 2 * jt + 1) /4. / pi)\
                                          * CG(2 * it, 0, 4, 0, 2 * jt, 0).doit()**2 / (2 * 2 * jt + 1) * Fct4(r)
                                
                                # projectile  
                                if  it == jt and it2 == jt2:
                                    if IVIBROTP == 0:
                                        if ip != jp:
                                            cc = betcahvp[ip][jp] * Fcp(r)
                                            if BetaP != 0.:
                                                cc = cc / BetaP
                                            
                                        else:
                                            cc = betcahvp[ip][jp] * Fcp2v(r) / sqrt(4. * pi)
                                        C = C + cc

                                    else:

                                        C = C + sqrt((2 * 2 * ip + 1) * 5 * (2 * 2 * jp + 1) /4. / pi)\
                                          * CG(2 * ip, 0, 2, 0, 2 * jp, 0).doit()**2 / (2 * 2 * jp + 1) * Fcp2(r)\
                                          * (betp[ip][jp] + 2. * sqrt(5./pi) * betp[ip][jp]**2/7.)

                                        C = C + sqrt((2 * 2 * ip + 1) * 9 * (2 * 2 * jp + 1) /4. / pi)\
                                          * CG(2 * ip, 0, 4, 0, 2 * jp, 0).doit()**2 / (2 * 2 * jp + 1) * Fcp4(r)
                                
                                if  ip == jp and it == jt:
                                    if it2 != jt2:
                                        cc = betcahv2[it2][jt2] * Fctt(r)
                                        if BetaT2 != 0.:
                                            cc = cc / BetaT2
                                    else:
                                        cc = betcahv2[it2][jt2] * Fct2v(r) / sqrt(4. * pi)

                                    C = C + cc
                                
                                # Excitation energy
                                if it == jt and ip == jp and it2 == jt2:
                                    if IVIBROTT == 0:
                                        C = C + omeahv[it]
                                    else:
                                        C + C + erott[it]
                                    
                                    if IVIBROTP == 0:
                                        C = C + omeahvp[ip]
                                    else:
                                        C = C + erotp[ip]
                                    C = C + omeahv2[it2]
                                
                                cpot[i][j] = C

        C0 = cpot[0][0]
        for i in range(0, Ndim):
            cpot[i][i] = cpot[i][i] - Vn(r)
            #cpot[i][i] = cpot[i][i] - C0
    # transfer coupling 
    if Ntrans == 1:
        cpot[1][Nlevel] = Ftrans(r)
        cpot[Nlevel][1] = Ftrans(r)
        cpot[Nlevel][Nlevel] = cpot[1][1] - Qtrans + (Zpro + iq) * (Ztar - iq) / r * Hbar / 137. - Zpro * Ztar / r * Hbar /137.
    return cpot

def Rkuta(psi0, phi0):
  """_summary_
  Wave functions at R=RMIN+DR
  """
  Ak1 = np.zeros((Nlevelmax), dtype=complex)
  Ak2 = np.zeros((Nlevelmax), dtype=complex)
  Ak3 = np.zeros((Nlevelmax), dtype=complex)
  Ak4 = np.zeros((Nlevelmax), dtype=complex)
  Bk1 = np.zeros((Nlevelmax), dtype=complex)
  Bk2 = np.zeros((Nlevelmax), dtype=complex)
  Bk3 = np.zeros((Nlevelmax), dtype=complex)
  Bk4 = np.zeros((Nlevelmax), dtype=complex)

  psi1 = np.zeros((Nlevelmax), dtype=complex)  
  Fac = dr * (2. * ReduceMass / Hbar**2)
  R_i  = R_min
  R_i1 = R_min + dr/2.
  R_i2 = R_min + dr

  for i0 in range(Nlevel):
      for ic in range(Nlevel):
          Ak1[i0] = Ak1[i0] + Fac*CPOT[i0][ic][0]*psi0[ic]
      Ak1[i0] = Ak1[i0] - Fac*(E_i - V(R_i, L))*psi0[i0]
      Bk1[i0] = dr * phi0[i0]

  for i0 in range(Nlevel):
      for ic in range(Nlevel):
          Ak2[i0] = Ak2[i0] + Fac*CPOTH[i0][ic] * (psi0[ic] + 1./2. *Bk1[ic])
      Ak2[i0] = Ak2[i0] - Fac*(E_i - V(R_i1, L)) * (psi0[i0] + 1./2. *Bk1[i0])
      Bk2[i0] = dr * (phi0[i0] + 1./2. * Ak1[i0])

  for i0 in range(Nlevel):
      for ic in range(Nlevel):
          Ak3[i0] = Ak3[i0] + Fac*CPOTH[i0][ic] * (psi0[ic] + 1./2. *Bk2[ic])
      Ak3[i0] = Ak3[i0] - Fac*(E_i - V(R_i1, L)) * (psi0[i0] + 1./2. *Bk2[i0])
      Bk3[i0] = dr * (phi0[i0] + 1./2. * Ak2[i0])
  
  for i0 in range(Nlevel):
      for ic in range(Nlevel):
        Ak4[i0] = Ak4[i0] + Fac*CPOT[i0][ic][1] * (psi0[ic] + Bk3[ic])
      Ak4[i0] = Ak4[i0] - Fac*(E_i - V(R_i2, L)) * (psi0[i0] + Bk3[i0])
      Bk4[i0] = dr * (phi0[i0] + Ak3[i0])

  for i0 in range(Nlevel):
      psi1[i0] = psi0[i0] + 1./6. * (Bk1[i0] + 2*Bk2[i0] + 2*Bk3[i0] + Bk4[i0])

  return psi1

def stabilize(xi1, xi, aa, ir1):
    # Initialize arrays
    psi  = np.zeros((Nlevel, Nlevel), dtype=complex)
    cc   = np.zeros((Nlevel, Nlevel), dtype=complex)
    cin  = np.zeros((Nlevel, Nlevel), dtype=complex)
    aa0  = np.zeros((Nlevel, Nlevel), dtype=complex)
    xid  = np.zeros((Nlevel, Nlevel), dtype=complex)
    xid1 = np.zeros((Nlevel, Nlevel), dtype=complex)

    # Copy matrices
    for i in range(Nlevel):
        for j in range(Nlevel):
            aa0[i, j] = aa[i, j]
            xid[i, j] = xi[i, j]
            xid1[i, j] = xi1[i, j]

    # Calculate factor
    fac = dr**2 * (2.0 * ReduceMass / Hbar**2)

    # Transform to physical wave functions psi from xi
    r = R_min + ir1 * dr
    for i0 in range(Nlevel):
        for ic in range(Nlevel):
            cc[i0, ic] = -fac / 12.0 * CPOT[i0, ic, ir1]
            if i0 == ic:
                cc[i0, ic] += -fac / 12.0 * (V(r, L) - E_i) + 1.0

    cin = np.linalg.inv(cc)  # Matrix inversion
    for ich in range(Nlevel):
        for i0 in range(Nlevel):
            psi[i0, ich] = 0.0
            for ic in range(Nlevel):
                psi[i0, ich] += cin[i0, ic] * xi[ic, ich]

    # Update matrix aa
    for i in range(Nlevel):
        for j in range(Nlevel):
            aa[i, j] = 0.0
            for k in range(Nlevel):
                aa[i, j] += psi[i, k] * aa0[k, j]

    # Update xi and xi1
    cin = np.linalg.inv(psi)  # Matrix inversion
    for i in range(Nlevel):
        for j in range(Nlevel):
            xi[i, j] = 0.0
            xi1[i, j] = 0.0
            for k in range(Nlevel):
                xi[i, j] = xi[i, j] + xid[i, k] * cin[k, j]
                xi1[i, j] = xi1[i, j] + xid1[i, k] * cin[k, j]

    return xi1, xi, aa

def Numerov():
    """_summary_
        Subroutine for integration of the c.c. eqs. by modified 
    """
    # Initialize arrays
    psi    = np.zeros(Nlevel, dtype=complex)
    psi0   = np.zeros(Nlevel, dtype=complex)
    psi1   = np.zeros(Nlevel, dtype=complex)
    phi0   = np.zeros(Nlevel, dtype=complex)
    xi     = np.zeros((Nlevel, Nlevel), dtype=complex)
    xi0    = np.zeros((Nlevel, Nlevel), dtype=complex)
    xi1    = np.zeros((Nlevel, Nlevel), dtype=complex)
    aa     = np.zeros((Nlevel, Nlevel), dtype=complex)
    bb     = np.zeros((Nlevel, Nlevel), dtype=complex)
    bin    = np.zeros((Nlevel, Nlevel), dtype=complex)
    bb0    = np.zeros((Nlevel, Nlevel), dtype=complex)
    bb2    = np.zeros((Nlevel, Nlevel), dtype=complex)
    bb20   = np.zeros((Nlevel, Nlevel), dtype=complex)
    cc     = np.zeros((Nlevel, Nlevel), dtype=complex)
    cin    = np.zeros((Nlevel, Nlevel), dtype=complex)
    dd0    = np.zeros((Nlevel, Nlevel), dtype=complex)
    dd1    = np.zeros((Nlevel, Nlevel), dtype=complex)
    dd     = np.zeros(Nlevel, dtype=complex)
    ech    = np.zeros(Nlevel, dtype=float)
    ech2   = np.zeros(Nlevel, dtype=float)
    fcw    = np.zeros(201)
    gcw    = np.zeros(201)
    fpcw   = np.zeros(201)
    gpcw   = np.zeros(201)
    sigmad = np.zeros(201)
    iexp   = np.zeros(201, dtype=int)

    # Constants
    ai       = 1j  # Imaginary unit
    fac      = dr**2 * (2.0 * ReduceMass / Hbar**2)
    iterat   = R_iterat#int((R_max - R_min) / dr)
    #R_max    = R_min + iterat * dr
    ibarrier = int((R_barrier - R_min) / dr + 1e-6)


    # Initial conditions
    for io in range(Nlevel):
        aa[io, io] = 1.0
        for j1 in range(Nlevel):
            psi[j1]  = 0.0
            psi0[j1] = 0.0
            psi1[j1] = 0.0
            phi0[j1] = 0.0

        if io == 0:
            for io2 in range(Nlevel):
                ech[io2]  = E_i - V(R_min, L) - CPOT[io2, io2, 0] 
                ech2[io2] = E_i - CPOT[io2, io2, -2] 
                
        if ech[io] > 0.0:
            k = np.sqrt(2.0 * ReduceMass  / Hbar**2 * ech[io])
            psi0[io] = np.exp(-ai * k * R_min)
            phi0[io] = -ai * k * psi0[io]
        else:
            k = np.sqrt(2.0 * ReduceMass  / Hbar**2 * abs(ech[io]))
            psi0[io] = np.exp(k * R_min)
            phi0[io] = k * psi0[io]

        # Call Runge-Kutta subroutine
        psi1 = Rkuta(psi0, phi0)

        for i0 in range(Nlevel):
            xi0[i0, io] = (1.0 - fac / 12.0 * (V(R_min, L) - E_i)) * psi0[i0]
            xi1[i0, io] = (1.0 - fac / 12.0 * (V(R_min + dr, L) - E_i)) * psi1[i0]
            for ic in range(Nlevel):
                xi0[i0, io] = xi0[i0, io] - fac / 12.0 * CPOT[i0, ic, 0] * psi0[ic]
                xi1[i0, io] = xi1[i0, io] - fac / 12.0 * CPOT[i0, ic, 1] * psi1[ic]

    # Iterations   iterat + 2
    for ir in range(2, iterat + 2):
        r  = R_min + dr * ir
        r0 = R_min + dr * (ir - 2)
        r1 = R_min + dr * (ir - 1)

        # Matrix calculations
        for i0 in range(Nlevel):
            for ic in range(Nlevel):
                dd0[i0, ic] = fac / np.sqrt(12.0) * CPOT[i0, ic, ir - 1]
                if i0 == ic:
                    dd0[i0, ic] = dd0[i0, ic] + fac / np.sqrt(12.0) * (V(r1, L) - E_i) + np.sqrt(3.0)

        for i0 in range(Nlevel):
            for ic in range(Nlevel):
                dd1[i0, ic] = -1.0 if i0 == ic else 0.0
                for ik in range(Nlevel):
                    dd1[i0, ic] = dd1[i0, ic] + dd0[i0, ik] * dd0[ik, ic]

        for ich in range(Nlevel):
            for i0 in range(Nlevel):
                xi[i0, ich] = -xi0[i0, ich]
                for ic in range(Nlevel):
                    xi[i0, ich] = xi[i0, ich] + dd1[i0, ic] * xi1[ic, ich]

        if ir == iterat + 1:
            break
        
        if ir == ibarrier:
            xi1, xi, aa = stabilize(xi1, xi, aa, ir)

        for ich in range(Nlevel):
            for i0 in range(Nlevel):
                xi0[i0, ich] = xi1[i0, ich]
                xi1[i0, ich] = xi[i0, ich]
    
    # Matching to Coulomb wave function at R_max
    for io in range(Nlevel):
        # Matrix calculations
        for i0 in range(Nlevel):
            for ic in range(Nlevel):
                cc[i0, ic] = -fac / 12.0 * CPOT[i0, ic, iterat - 1]
                if i0 == ic:
                    cc[i0, ic] = cc[i0, ic] - fac / 12.0 * (V(R_max - dr, L) - E_i) + 1.0
                
        cin = np.linalg.inv(cc)

        for i0 in range(Nlevel):
            psi0[i0] = 0.0
            for ic in range(Nlevel):
                psi0[i0] = psi0[i0] + cin[i0, ic] * xi0[ic, io]

        for i0 in range(Nlevel):
            for ic in range(Nlevel):
                cc[i0, ic] = -fac / 12.0 * CPOT[i0, ic, iterat + 1]
                if i0 == ic:
                    cc[i0, ic] = cc[i0, ic] - fac / 12.0 * (V(R_max + dr, L) - E_i) + 1.0
                

        cin = np.linalg.inv(cc)

        for i0 in range(Nlevel):
            psi[i0] = 0.0
            for ic in range(Nlevel):
                psi[i0] = psi[i0] + cin[i0, ic] * xi[ic, io]

        # Additional calculations for Coulomb wave function matching
        # (Omitted for brevity)
        for ii in range(Nlevel):
            # Calculate energy difference
            ec = E_i - CPOT[ii, ii, iterat-1]
            
            if ec < 0.0:
                # For negative energy (exponential solutions)
                r1 = R_max - dr
                r2 = R_max + dr
                ak = np.sqrt(2.0 * ReduceMass * abs(ec) / Hbar**2)

                # Calculate bb0 and bb20
                bb0[ii, io]  =  (np.exp(-ak * r2) * psi0[ii] - np.exp(-ak * r1) * psi[ii]) / (np.exp(ak * (r1 - r2))) - np.exp(-ak * (r1 - r2))
                bb20[ii, io] = -(np.exp(ak * r2) * psi0[ii] - np.exp(ak * r1) * psi[ii]) / (np.exp(ak * (r1 - r2)) - np.exp(-ak * (r1 - r2)))
            else:
                # For positive energy (Coulomb wave functions)
                rho = np.sqrt(2.0 * ReduceMass * ec) / Hbar * (R_max - dr)
                eta = (Zpro * Ztar / 137.0) * np.sqrt(ReduceMass / (2.0 * ec))

                fcw[L] = coulombf(L, eta, rho)
                gcw[L] = coulombg(L, eta, rho)

                # Calculate Coulomb wave functions
                cwup0 = gcw[L] + ai * fcw[L]
                cwdown0 = gcw[L] - ai * fcw[L]

                # Recalculate for rmax + dr
                ec = E_i - CPOT[ii, ii, -1]  # cpot(ii,ii,iterat+1)
                rho = np.sqrt(2.0 * ReduceMass * ec) / Hbar * (R_max + dr)
                eta = (Zpro * Ztar / 137.0) * np.sqrt(ReduceMass / (2.0 * ec))

                fcw[L] = coulombf(L, eta, rho)
                gcw[L] = coulombg(L, eta, rho)
                cwup1 = gcw[L] + ai * fcw[L]
                cwdown1 = gcw[L] - ai * fcw[L]

                # Calculate bb0 and bb20
                bb0[ii, io] = (cwup0 * psi[ii] - cwup1 * psi0[ii]) / (cwup0 * cwdown1 - cwup1 * cwdown0)
                bb20[ii, io] = (cwdown1 * psi0[ii] - cwdown0 * psi[ii]) / (cwup0 * cwdown1 - cwup1 * cwdown0)

    # Calculate bb and bb2 matrices
    for i in range(Nlevel):
        for j in range(Nlevel):
            bb[i, j] = 0.0
            bb2[i, j] = 0.0
            for j2 in range(Nlevel):
                bb[i, j] = bb[i, j] + bb0[i, j2] * aa[j2, j]
                bb2[i, j]= bb2[i, j] + bb20[i, j2] * aa[j2, j]

    # Calculate the inverse of bb
    bin = np.linalg.inv(bb)
    # Penetration probability calculation
    P = 0.0
    ref = 0.0
    for io in range(Nlevel):
        if ech[io] < 0.0:
            break
        k  = np.sqrt(2.0 * ReduceMass / Hbar**2 * ech[io])
        kk = np.sqrt(2.0 * ReduceMass / Hbar**2 * E_i)
        P  = P + abs(bin[io, 0])**2 * k / kk

    for io in range(Nlevel):
        if ech2[io] < 0.0:
            continue
        k2 = np.sqrt(2.0 * ReduceMass / Hbar**2 * ech2[io])
        kk = np.sqrt(2.0 * ReduceMass / Hbar**2 * E_i)
        dummy2 = 0.0
        for io2 in range(Nlevel):
            dummy2 = dummy2 + bb2[io, io2] * bin[io2, 0]
        ref = ref + abs(dummy2)**2 * k2 / kk

    #if P > 1.0:
    #  print('Penetrability > 1 !!!!!')
    #  print('Something must be strange :^(')

    return P


"""
Principal physical constants

"""
Hbar = 197.329E0
Pi = 3.141592
Amu = 938E0
e2 = 1.44
Element = ['O','H','He','Li','Be','B ','C','N','O','F','Ne',
          'Na','Mg','Al','Si','P','S','Cl','Ar','K ','Ca','Sc','Ti','V ',
          'Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
          'Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In',
          'Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm',
          'Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re',
          'Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra',
          'Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md',
          'No','Lr','XX','X1','X2','X3','X4','04']


"""Principal global variables
P                   - penetrability 
Sigma               - fusion cross section, unit mb
Spin                - mean angular momentum
  
Apro, Zpro, Rpro    - atomic #, proton # and radius of the projectile
Atar, Ztar, Rtar    - those of the target
ReduceMass          - reduced mass, unit MeV/c**2
Amu                 - nucleon mass, unit MeV/c**2
E                   - bombarding energy in the center of mass frame, unit MeV
  
V0,R0,A0            - depth, range, and dissuseness parameters of uncoupled 
                      nuclear potential, which is assumed to be a Woods-Saxon form
  
IVIBROTT (IVIBROTP) - option for intrinsic degree of freedom 
                      = -1; no excitation (inert)
                      =  0; vibrational coupling 
                      =  1; rotational coupling

Ntar (Npro)         - the number of levels to be included    

BetaTn (BetaPn)     - i am not very clear about it (Liao)
BetaT (BetaP)       - defeormation parameter
OmegaT (OmegaP)     - excitation energy of the oscillator
LambdaT (LambdaP)   - multipolarity
NphononT (NphononP) - the number of phonon to be included
Beta2T (Beta2P)     - static quadrupole deformation parameter
Beta4T (Beta4P)     - static hexadecapole deformation parameter
BetaT2n             - is the same as BetaTn but for the second mode
BetaT2              - is the same as BETAT but for the second mode
OmegaT2             - is the same as OmegaT but for the second mode
LambdaT2            - is the same as LambdaT but for the second mode
NphononT2           - is the same as NphononT but for the second mode

E2T (E2P)           - excitation energy of 2+ state in a rotational band
NrotT (NrotP)       - the number of levels in the rotational band to be included 
                      (up to I^pi=2*NROT+ states are included)

L                   - angular momentum of the relative motion
CPOT                - coupling matrix
"""
start_time = time.time()

# collective motion coupled parameter
with open("ccfull.inp", "r") as f:
    values = [
        float(x.strip()) if "." in x else int(x.strip())
        for line in f if line.strip()
        for x in line.split(",")
    ]

Apro, Zpro, Atar, Ztar       = values[0], values[1], values[2], values[3]
R0P, IVIBROTP, R0T, IVIBROTT = values[4], values[5], values[6], values[7]

Rpro = R0P*Apro**(1./3)
Rtar = R0T*Atar**(1./3)
ReduceMass = Apro*Atar/(Apro+Atar) * Amu
h2m = Hbar ** 2 / (2 * ReduceMass)

Ntar, Npro = 0, 0
# Mode of excitation for target
OmegaT, BetaT, LambdaT, NphononT      = 0., 0., 0., 0
OmegaT2, BetaT2, LambdaT2, NphononT2  = 0., 0., 0., 0
E2T, Beta2T, Beta4T, NrotT            = 0., 0., 0., 0
if IVIBROTT == -1:  # for inert
    Ntar = 0
elif IVIBROTT == 0:  # for vibration
    OmegaT, BetaT, LambdaT, NphononT  = values[8], values[9], values[10], values[11]
    Ntar = NphononT
    # Second excitation for target
    OmegaT2, BetaT2, LambdaT2, NphononT2  = values[12], values[13], values[14], values[15]
elif IVIBROTT == 1:  # for rotation
    E2T, Beta2T, Beta4T, NrotT        = values[8], values[9], values[10], values[11]
    Ntar = NrotT
else:
    print(f"Invalid input for IVIBROTT: {IVIBROTT}. Expected values are -1, 0, or 1.")
    exit()


# Mode of excitation for projectile
OmegaP, BetaP, LambdaP, NphononP  = 0., 0., 0., 0.
E2P, Beta2P, Beta4P, NrotP        = 0., 0., 0., 0.
if IVIBROTP == -1:  # for inert
    Npro = 0
elif IVIBROTP == 0:  # for vibration
    OmegaP, BetaP, LambdaP, NphononP  = values[16], values[17], values[18], values[19]
    Npro = NphononP
elif IVIBROTP == 1:  # for rotation
    E2P, Beta2P, Beta4P, NrotP        = values[16], values[17], values[18], values[19]
    Npro = NrotP
else:
    print(f"Invalid input for IVIBROTP: {IVIBROTP}. Expected values are -1, 0, or 1.")
    exit()

# transfer coupled parameter
Ntrans, Qtrans, Ftr =  values[20], values[21], values[22]
iq = 0
# Potential parameters
V0, R0, A0 = values[23], values[24], values[25]

# Energy interval
Emin, Emax, dE = values[26], values[27], values[28]
E = np.arange(Emin, Emax, dE)
# Radial interval
dr = 0.0001
R = np.arange(1, 30, dr)

# beta^n_lambda and beta^C_lambda
BetaTn, BetaPn = 1., 1.
BetaT2n, BetaP2n = 1., 1.


# initialzation
Nlevelmax = 30
Nlevel = (Npro + 1) * (Ntar + 1) + Ntrans
if Nlevel > Nlevelmax:
    print(f" too many channels, Nlevels = {Nlevel} is lager than Nlevelmax {Nlevelmax}")
    exit()
    
imutual  = [[0 for _ in range(Nlevelmax + 1)] for _ in range(Nlevelmax + 1)]
Ev       = np.zeros((Nlevelmax), dtype=float)
Evec     = np.zeros((Nlevelmax, Nlevelmax), dtype=float)

betnahv  = [[0.0 for _ in range(Nlevelmax + 1)] for _ in range(Nlevelmax + 1)]
betcahv  = [[0.0 for _ in range(Nlevelmax + 1)] for _ in range(Nlevelmax + 1)]
omeahv   = [0.0 for _ in range(Nlevelmax + 1)]
betnahv2 = [[0.0 for _ in range(Nlevelmax + 1)] for _ in range(Nlevelmax + 1)]
betcahv2 = [[0.0 for _ in range(Nlevelmax + 1)] for _ in range(Nlevelmax + 1)]
omeahv2  = [0.0 for _ in range(Nlevelmax + 1)]
betnahvp = [[0.0 for _ in range(Nlevelmax + 1)] for _ in range(Nlevelmax + 1)]
betcahvp = [[0.0 for _ in range(Nlevelmax + 1)] for _ in range(Nlevelmax + 1)]
omeahvp  = [0.0 for _ in range(Nlevelmax + 1)]

bett     = [[0.0 for _ in range(Nlevelmax + 1)] for _ in range(Nlevelmax + 1)]
erott    = [0.0 for _ in range(Nlevelmax + 1)]
betp     = [[0.0 for _ in range(Nlevelmax + 1)] for _ in range(Nlevelmax + 1)]
erotp    = [0.0 for _ in range(Nlevelmax + 1)]

# Show and record important parameters
SystemName = str(Apro) + str(Element[Zpro]) + '+' + str(Atar) + str(Element[Ztar])
Filename_output = "CC_python_output_"+SystemName +".dat"
Filename_sigma  = "CC_python_sigma_"+SystemName +".dat"
print(f"Reaction system: {Apro}{Element[Zpro]} + {Atar}{Element[Ztar]}")
print(f"Simulation range E from {E[0]} MeV To {E[-1]} MeV, dE = {dE} MeV")
print(f"Simulation range R from {R[0]:.2f} fm To {R[-1]:.2f} fm, dR = {dr} fm")
print(f"------------------------------")
Pot_para1, Pot_para2, Pot_para3 = WS_parameter()
print(f"Potential parameters: V0=  {Pot_para1:.2f}(MeV), a= {Pot_para2:.2f}(fm), r0= {Pot_para3:.2f}(fm), R= {Pot_para3*(Apro**(1./3)+Atar**(1./3)):.2f}(fm)")
print(f"------------------------------")
print(f"Reaction system: {Apro}{Element[Zpro]} + {Atar}{Element[Ztar]}", file=open(Filename_output, "w"))
print(f"Simulation range E from {E[0]} MeV To {E[-1]} MeV, dE = {dE} MeV", file=open(Filename_output, "a"))
print(f"Simulation range R from {R[0]:.2f} fm To {R[-1]:.2f} fm, dR = {dr} fm", file=open(Filename_output, "a"))
print(f"Potential parameters: V0=  {Pot_para1:.2f}(MeV), a= {Pot_para2:.2f}(fm), r= {Pot_para3:.2f}(fm), r0= {Pot_para3*(Apro**(1./3)+Atar**(1./3)):.2f}(fm)", file=open(Filename_output, "a"))
print(f"------------------------------")

#Coulomb barrier
local_max_indices = argrelmax(V(R, 0))[0]
local_min_indices = argrelmin(V(R, 0))[0]
V_barrier = V(R, 0)[local_max_indices].item()
R_barrier = R[local_max_indices].item()
V_bottom = V(R, 0)[local_min_indices].item()
R_bottom = R[local_min_indices].item()
print(f"Coulomb barrier idex: {local_max_indices}")
print(f"Coulomb barrier postion: {R_barrier:.4f} fm")
print(f"Coulomb barrier energy: {V_barrier:.4f} MeV")
print(f"Coulomb bottom idex: {local_min_indices}")
print(f"Coulomb bottom position: {R_bottom:.4f} fm")
print(f"Coulomb bottom energy: {V_bottom:.4f} MeV")

print(f"------------------------------",            file=open(Filename_output, "a"))
print(f"Coulomb barrier idex: {local_max_indices}", file=open(Filename_output, "a"))
print(f"Coulomb barrier postion: {R_barrier:.4f} fm",   file=open(Filename_output, "a"))
print(f"Coulomb barrier energy: {V_barrier:.4f} MeV",   file=open(Filename_output, "a"))
print(f"Coulomb bottom idex: {local_min_indices}",  file=open(Filename_output, "a"))
print(f"Coulomb bottom position: {R_bottom:.4f} fm",    file=open(Filename_output, "a"))
print(f"Coulomb bottom energy: {V_bottom:.4f} MeV",     file=open(Filename_output, "a"))
print(f"----------------------------",              file=open(Filename_output, "a"))


print(f"------------------------------")
print(f"Mode of excitation for target")
print(f"Mode of excitation for target", file=open(Filename_output, "a"))
if Ntar != 0:
    if IVIBROTT == 0:
        # Phonon excitation in the target
        print(f"Phonon Excitation in the targ.: beta={BetaT:.3f}, omega={OmegaT:.2f} (MeV), Lambda={LambdaT}, Nph={NphononT}")
        BetaTn = BetaT
        print("    ")
        ans = input("Different beta_N from beta_C for this mode (n/y)? ").strip().lower()
        if ans == 'y':
            BetaTn = float(input("beta_N=? "))
        
        # Write to file (assuming file handle `f` is opened elsewhere)
        print(f"Phonon Excitation in the targ.: beta_N={BetaTn:.3f}, beta_C={BetaT:.3f}, r0={R0T:.2f} (fm),\n"
              f"                              omega={OmegaT:.2f} (MeV), Lambda={LambdaT}, Nph={NphononT}", file=open(Filename_output, "a"))

    if IVIBROTT == 1:
        # Rotational excitation in the target
        output_str = (f"Rotational Excitation in the targ.: beta2={Beta2T:.3f}, beta4={Beta4T:.3f}, r0={R0T:.2f} (fm),\n"
                      f"                                   E2={E2T:.2f} (MeV), Nrot={NrotT}")

        print(output_str, file=open(Filename_output, "a"))

if NphononT2 != 0:
    # Phonon excitation in the target (second mode)
    print(f"Phonon Excitation in the targ.: beta={BetaT2:.3f}, omega={OmegaT2:.2f} (MeV), Lambda={LambdaT2}, Nph={NphononT2}")
    
    BetaT2n = BetaT2
    print("    ")
    ans = input("Different beta_N from beta_C for this mode (n/y)? ").strip().lower()
    if ans == 'y':
        BetaT2n = float(input("beta_N=? "))
    
    # Write to file (assuming file handle `f` is opened elsewhere)
    print(f"Phonon Excitation in the targ.: beta_N={BetaT2n:.3f}, beta_C={BetaT2:.3f}, r0={R0T:.2f} (fm),\n"
          f"                              omega={OmegaT2:.2f} (MeV), Lambda={LambdaT2}, Nph={NphononT2}", file=open(Filename_output, "a"))
    mutual()

print(f"------------------------------")
print(f"Mode of excitation for projectile")
print(f"Mode of excitation for projectile", file=open(Filename_output, "a"))

if Npro != 0:
    if IVIBROTP == 0:
        # Phonon excitation in the projectile
        print(f"Phonon Excitation in the proj.: beta={BetaP:.3f}, omega={OmegaP:.2f} (MeV), Lambda={LambdaP}, Nph={NphononP}")
        BetaPn = BetaP
        print("    ")
        ans = input("Different beta_N from beta_C for this mode (n/y)? ").strip().lower()
        if ans == 'y':
            betapn = float(input("beta_N=? "))
        
        # Write to file (assuming file handle `f` is opened elsewhere)
        print(f"Phonon Excitation in the proj.: beta_N={BetaPn:.3f}, beta_C={BetaP:.3f}, r0={R0P:.2f} (fm),\n"
              f"                              omega={OmegaP:.2f} (MeV), Lambda={LambdaP}, Nph={NphononP}", file=open(Filename_output, "a"))

    if IVIBROTP == 1:
        # Rotational excitation in the projectile
        output_str = (f"Rotational Excitation in the proj.: beta2={Beta2P:.3f}, beta4={Beta4P:.3f}, r0={R0P:.2f} (fm),\n"
                      f"                                   E2={E2P:.2f} (MeV), Nrot={NrotP}")
        print(output_str)
        print(output_str, file=open(Filename_output, "a"))

print(f"------------------------------")
print(f"Transfer coupled mode")
print(f"Transfer coupled mode", file=open(Filename_output, "a"))
if Ntrans != 0:
    print(f"Transfer channel: Strength= {Ftr}, Q = {Qtrans} MeV")
    print(f"Transfer channel: Strength= {Ftr}, Q = {Qtrans} MeV", file=open(Filename_output, "a"))

grotation()
anharmonicity()


df = pd.DataFrame(bett)
print('bett:', file=open(Filename_output, "a"))
print(df.iloc[:6, :6], file=open(Filename_output, "a"))

df = pd.DataFrame(betcahvp)
print('betcahvp:', file=open(Filename_output, "a"))
print(df.iloc[:6, :6], file=open(Filename_output, "a"))


print(f'Nlevel = {Nlevel}')
print(f'Nlevel = {Nlevel}', file=open(Filename_output, "a"))
dr = 0.05
R_iterat = int((R[-1] - R_bottom)/dr)
R_min, R_max = R_bottom, R_bottom + dr*R_iterat

print(f"Rmin = {R_min:.4f} fm, Rmax = {R_max:.4f} fm, iterat times = {R_iterat}")
print(f"Rmin = {R_min:.4f} fm, Rmax = {R_max:.4f} fm, iterat times = {R_iterat}", file=open(Filename_output, "a"))
print('---------------------------------')
print('Initialize the Hamiltonian matrix0')
coupled_matrix0()
CPOT  = np.zeros((Nlevelmax, Nlevelmax, R_iterat+2), dtype=float)
CPOT0 = np.zeros((Nlevelmax, Nlevelmax), dtype=float)
CPOTH = np.zeros((Nlevelmax, Nlevelmax), dtype=float)
print('---------------------------------')
print('Initialize the Hamiltonian matrix at R')
for ir in range(R_iterat + 2):
    CPOT0 = CoupledMatrix(R_min + dr*ir)
    
    df = pd.DataFrame(CPOT0)
    print(f'ir = {ir} and r = {R_min + dr*ir:.5f}fm', file=open(Filename_output, "a"))
    print(f'ir = {ir} and r = {R_min + dr*ir}fm', file=open(Filename_output, "a"))
    print('coupled_matrix:', file=open(Filename_output, "a"))
    print(df.iloc[:Nlevel, :Nlevel], file=open(Filename_output, "a"))

    for i in range(Nlevel):
        for j in range(Nlevel):
            CPOT[i][j][ir] = CPOT0[i][j]


print('---------------------------------')
R_h = R_min + dr/2.
CPOTH = CoupledMatrix(R_h)
print(f'Initialize the Hamiltonian matrix at R+dr/2, rh = {R_h:.4f}fm')

print('---------------------------------')
print('Begin to calculate the penetration probability at certain E and L')
CalData = []

# Energy loop
for ie in range(len(E)):
    E_i = E[0] + ie * dE
    
    sigma = 0.
    spin = 0.
    P0 = 0.
    # Angular momentum loop
    for L in range(0, 201):
        
        # find the V(E, L) coulomb barrier postion
        local_max_indices = argrelmax(V(R, L))[0]
        local_min_indices = argrelmin(V(R, L))[0]
        if V(R, L)[local_min_indices].item() > E_i or R[local_min_indices].item() < 0:
            P = 0.
            #print(f" can not find a barrier when  E = {E_i} MeV and L > {L} hbar")
            break
        
        P = Numerov()
        if L == 0:
            P0 = P
        sigma = sigma + (2.*L+1)*P*pi*Hbar**2/2.0/ReduceMass/E_i*10.
        spin  = spin  + (2.*L+1)*P*pi*Hbar**2/2.0/ReduceMass/E_i*10.*L
        if (2.*L+1)*P*pi*Hbar**2/2.0/ReduceMass/E_i*10. < sigma * 0.0001:
            break
    if sigma < 0.01:
      print(f' E = {E_i:1f} MeV,   Sigma = {sigma:.5e} mb, <L> = {spin/sigma:.5f} hbar, <P0> = {P0:.5e}')
    else:
      print(f' E = {E_i:1f} MeV,   Sigma = {sigma:.5f} mb, <L> = {spin/sigma:.5f} hbar, <P0> = {P0:.5f}')
    CalData.append([E_i, sigma, spin, P0])



# Output the fusion exciation function
with open(Filename_sigma, "w") as file:
    # Insert header
    print("E_i (MeV),    Sigma (mb),     <L> (hbar),       <P0>", file=file)
    for E_i, sigma, spin, P0 in CalData:
        if sigma < 0.01:
          print(f"{E_i}, {sigma:.5e}, {spin/sigma:.5f}, {sigma:.5e}", file=file)
        else:
          print(f"{E_i}, {sigma:.5f}, {spin/sigma:.5f}, {sigma:.5f}", file=file)

print("Data has been output into sigma.dat")
end_time = time.time()
time_difference = end_time - start_time
print(f"The total time taken by the procedure: {time_difference} s / {time_difference/3600:.2f} hour")