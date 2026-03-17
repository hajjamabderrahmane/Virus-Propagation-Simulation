import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp, trapezoid, simpson , cumulative_trapezoid , odeint
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

# Chargement des données
data = {
    "Date": ["1/3/2020", "8/3/2020", "15/3/2020", "22/3/2020", "29/3/2020", "5/4/2020", 
             "12/4/2020", "19/4/2020", "26/4/2020", "3/5/2020", "10/5/2020", "17/5/2020", 
             "24/5/2020", "31/5/2020", "7/6/2020", "14/6/2020", "21/6/2020", "28/6/2020", 
             "5/7/2020", "12/7/2020", "19/7/2020", "26/7/2020", "2/8/2020", "9/8/2020", 
             "16/8/2020", "23/8/2020", "30/8/2020", "6/9/2020", "13/9/2020", "20/9/2020", 
             "27/9/2020", "4/10/2020", "11/10/2020", "18/10/2020", "25/10/2020", "1/11/2020", 
             "8/11/2020", "15/11/2020", "22/11/2020", "29/11/2020", "6/12/2020", "13/12/2020", 
             "20/12/2020", "27/12/2020", "3/1/2021", "10/1/2021", "17/1/2021", "24/1/2021", 
             "31/1/2021", "7/2/2021", "14/2/2021", "21/2/2021", "28/2/2021"],
    "I": [0, 21, 1034, 2580, 3560, 3782, 2887, 2342, 1003, 787, 457, 360, 
          467, 230, 197, 180, 172, 201, 155, 223, 175, 140, 135, 133, 140, 
          150, 186, 207, 263, 206, 382, 335, 281, 412, 362, 501, 788, 896, 
          1148, 1891, 1753, 2177, 2346, 2888, 3464, 3197, 2919, 2927, 2894, 
          1645, 1898, 2286, 2037],
    "S": [8250000, 8249976, 8248142, 8245613, 8244389, 8243842, 8245177, 
          8245852, 8248301, 8248642, 8249191, 8249360, 8249152, 8249600, 
          8249660, 8249687, 8249700, 8249630, 8249719, 8249580, 8249688, 
          8249745, 8249763, 8249760, 8249741, 8249723, 8249648, 8249617, 
          8249491, 8249298, 8249273, 8249374, 8249487, 8249216, 8249320, 
          8249053, 8248459, 8248300, 8247820, 8246981, 8246667, 8245873, 
          8245554, 8244465, 8243312, 8237916, 8244500, 8244480, 8244515, 
          8247020, 8246447, 8245702, 8246170],
    "R": [0, 3, 820, 1755, 1695, 1618, 1182, 1312, 396, 361, 238, 191, 
          327, 128, 114, 111, 106, 157, 117, 190, 127, 107, 100, 104, 114, 
          124, 160, 170, 243, 491, 341, 284, 225, 368, 308, 434, 744, 790, 
          1020, 1113, 1558, 1918, 2066, 2599, 3181, 8831, 2522, 2519, 2515, 
          1255, 1572, 1950, 1721],
    "D": [0, 0, 4, 52, 356, 758, 754, 494, 300, 210, 114, 89, 54, 42, 29, 
          22, 22, 12, 9, 7, 10, 8, 2, 3, 5, 3, 6, 6, 3, 5, 4, 7, 7, 4, 10, 
          12, 9, 14, 12, 15, 22, 32, 34, 48, 43, 56, 59, 74, 76, 80, 83, 
          62, 72]
}

df = pd.DataFrame(data)
df['Date'] = pd.to_datetime(df['Date'], dayfirst=True)
df['t'] = (df['Date'] - df['Date'].iloc[0]).dt.days

# Paramètres du modèle SIRD
r = 0.00000007927  # taux de contagion
a = 0.2172       # taux de guérison
b = 0.0157       # taux de mortalité

print("Estimation du nombre total de personnes infectées sur la période :")
df['dI_dt'] = r * df['S'] * df['I'] - (a + b) * df['I']
df['jours'] = (df['Date'] - df['Date'].iloc[0]).dt.days
# Vérification qu'il n'y a pas de NaN
print("Vérification des valeurs NaN :")
print(df.isna().sum())
# Intégration numérique
I_total_trapeze = trapezoid(df['dI_dt'], df['jours'])
I_total_simpson = simpson(df['dI_dt'], df['jours'])
# Résultats
print("\nRésultats d'intégration :")
print(f"Nombre total d'infections estimé (trapèzes): {I_total_trapeze:.0f}")
print(f"Nombre total d'infections estimé (Simpson): {I_total_simpson:.0f}")

print("Calcul du nombre cumulé de cas infectés et de guérisons :")
df['jours'] = (df['Date'] - df['Date'].iloc[0]).dt.days
# Calcul des intégrales cumulées
df['Nouveaux_cas'] = r * df['S'] * df['I']
df['Nouvelles_guerisons'] = a * df['I']
# Intégration numérique avec la méthode des trapèzes cumulés
df['I_cumule'] = cumulative_trapezoid(df['Nouveaux_cas'], df['jours'], initial=0)
df['R_cumule'] = cumulative_trapezoid(df['Nouvelles_guerisons'], df['jours'], initial=0)
# Résultats finaux
total_infections = df['I_cumule'].iloc[-1]
total_guerisons = df['R_cumule'].iloc[-1]
print(f"Nombre total d'infections cumulées: {total_infections:.0f}")
print(f"Nombre total de guérisons cumulées: {total_guerisons:.0f}")

print("Détermination de la durée totale de l'épidémie :")
# On cherche quand I(t) devient négligeable (par exemple < 1% du pic)
peak_I = df['I'].max()
threshold = 0.01 * peak_I
# Trouver le premier moment où I descend en dessous du seuil après le pic
peak_idx = df['I'].idxmax()
post_peak_df = df.loc[peak_idx:]
below_threshold = post_peak_df[post_peak_df['I'] < threshold]
if not below_threshold.empty:
    end_epidemic_idx = below_threshold.index[0]
    end_epidemic_date = df.loc[end_epidemic_idx, 'Date']
    duration = (end_epidemic_date - df['Date'].iloc[0]).days
    print(f"Durée estimée de l'épidémie: {duration} jours (jusqu'au {end_epidemic_date.date()})")
else:
    print("L'épidémie ne descend pas en dessous du seuil (I approche de 0 cas) dans les données disponibles alos il ya pas de fin d'épidemie en année étudiée . ")

print("Évaluation de l'impact des mesures de contrôle (simulation) :")
df['Date'] = pd.to_datetime(df['Date'], format='%d/%m/%Y')
start_date = df['Date'].min()
df['Days'] = (df['Date'] - start_date).dt.days
# --- Fonction pour calculer l'intégrale (méthode des trapèzes) ---
def calculate_integral(days, values):
    integral = 0
    for i in range(len(days)-1):
        delta = days[i+1] - days[i]
        avg_value = (values[i] + values[i+1]) / 2
        integral += avg_value * delta
    return integral


integral_avec_mesures = calculate_integral(df['Days'].values, df['I'].values)
# --- Modèle SIRD pour le scénario SANS mesures ---
N = 8_250_000  # Population totale (identique aux données)
# Paramètres SIRD (à ajuster pour correspondre à la dynamique initiale)
beta = 0.00000007927   # Taux de transmission
gamma = 0.2172   # Taux de guérison
mu = 0.0157     # Taux de mortalité

def SIRD_model(y, t, beta, gamma, mu):
    S, I, R, D = y
    dSdt = -beta * S * I / N
    dIdt = beta * S * I / N - gamma * I - mu * I
    dRdt = gamma * I
    dDdt = mu * I
    return [dSdt, dIdt, dRdt, dDdt]

# Conditions initiales (I=21 à t=7 jours, R=3, D=0)
I0 = 21
R0 = 3
D0 = 0
S0 = N - I0 - R0 - D0
y0 = [S0, I0, R0, D0]
# Temps (jours)
t = df['Days'].values
# Résolution du modèle SIRD
solution = odeint(SIRD_model, y0, t, args=(beta, gamma, mu))
S, I_sans_mesures, R, D = solution.T

# --- Calcul de l'intégrale SANS mesures ---
integral_sans_mesures = calculate_integral(t, I_sans_mesures)
# --- Impact des mesures ---
impact = integral_sans_mesures - integral_avec_mesures
# --- Résultats ---
print(f"Intégrale avec mesures: {integral_avec_mesures:.2f} cas-jours")
print(f"Intégrale sans mesures (SIRD): {integral_sans_mesures:.2f} cas-jours")
print(f"Impact des mesures: {impact:.2f} cas-jours évités")

print("Estimation du nombre de doses de vaccin nécessaires pour l'immunité collective :")
# Seuil d'immunité collective
S_c = (a + b) / r
P_c = 1 - S_c / y0[0]  # Proportion à vacciner
vaccines_needed = P_c * y0[0]

print(f"Nombre de doses de vaccin nécessaires pour l'immunité collective: {vaccines_needed:.0f}")

