# 🦠 Simulation de la propagation du COVID-19 – Modèle SIRD

## 📌 Description du projet

Ce projet a été réalisé dans le cadre du module **Analyse Numérique** à l’ENSAM Casablanca.
Il vise à modéliser et simuler la propagation du virus **COVID-19** à l’aide de méthodes numériques.

Le modèle utilisé est le **modèle SIRD**, qui permet de représenter l’évolution d’une épidémie dans une population.

---

## 🌍 Contexte : COVID-19

Le **COVID-19** est une maladie infectieuse causée par le coronavirus SARS-CoV-2, apparue en 2019.
Sa propagation rapide à l’échelle mondiale a nécessité l’utilisation de modèles mathématiques pour :

* Comprendre la dynamique de transmission
* Prédire l’évolution de l’épidémie
* Évaluer l’impact des mesures sanitaires

---

## ⚙️ Modélisation mathématique

Le modèle SIRD est défini par le système d’équations différentielles suivant :

* dS/dt = -r × S × I
* dI/dt = r × S × I - (a + b) × I
* dR/dt = a × I
* dD/dt = b × I

### 🔍 Interprétation :

* **S (Susceptible)** : individus sains
* **I (Infectés)** : individus contaminés par le COVID-19
* **R (Recovered)** : individus guéris et immunisés
* **D (Dead)** : individus décédés

---

## 🎯 Objectifs

* Modéliser la propagation du COVID-19 à partir de données réelles
* Appliquer les méthodes d’analyse numérique
* Étudier les indicateurs clés de l’épidémie :

  * Pic épidémique
  * Nombre de reproduction (R_0)
  * Seuil d’immunité collective

---

## 🧠 Méthodes utilisées

* Interpolation des données réelles
* Résolution d’équations non linéaires
* Dérivation numérique
* Intégration numérique
* Résolution d’équations différentielles

---

## 📊 Résultats attendus

* Courbes S(t), I(t), R(t), D(t)
* Analyse de la dynamique de propagation du COVID-19
* Comparaison entre modèle et données réelles

---

## 💻 Technologies utilisées

* Python
* NumPy
* Pandas
* Matplotlib

---

## 📌 Remarque

Ce projet est réalisé à des fins pédagogiques dans le cadre de l’apprentissage de l’analyse numérique et de la modélisation des épidémies.

---

## 📜 Licence

Projet académique – utilisation éducative uniquement.
