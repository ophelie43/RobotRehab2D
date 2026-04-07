
import pandas
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
import matplotlib.ticker as ticker # <--- Import important
from scipy.spatial import KDTree

file_path_cercle = "Data _position_CSV /CERCLE.csv" 
file_path_carre = "Data _position_CSV /Carre.csv"
file_path_aleatoire = "Data _position_CSV /suivi_robot.csv"


def extract_data_from_csv(file_path):
    df = pandas.read_csv(file_path)
    return df

def show_trajectories(df):

    # Extract data
    X_encodeur = df.X_Encodeur
    Y_encodeur = df.Y_Encodeur
    X_voulu = df.X_Voulu
    Y_voulu = df.Y_Voulu

    # Plot results
    fig, ax = plt.subplots()
    ax.grid()
    ax.plot(X_voulu, Y_voulu, color = 'green', label = 'Commande envoyée')
    ax.plot(X_encodeur, Y_encodeur, color = 'blue', label = 'Position réelle')
    ax.set_title("Commande envoyée vs position obtenue via lecture d'encodeur")
    ax.set_xlabel("Position en X [m]")
    ax.set_ylabel("Position en Y [m]")
    ax.legend()
    plt.show()

def compute_distance(reference: np.ndarray, prediction: np.ndarray):
    """
    Computes distance

    Parameters
    -----------
    reference: np.ndarray
        Reference value
    prediction : np.ndarray
        Calculated/predicted value
    Returns
    ---------
    distance : np.ndarray
    """
    ref = reference.T if reference.shape[0] == 2 else reference
    pred = prediction.T if prediction.shape[0] == 2 else prediction
    tree = KDTree(ref)
    distances, indices = tree.query(pred)

    return distances
def show_trajectories_with_distance(file_path):

    name1 = file_path.split("/")[1]
    name = name1.split(".")[0]
    df = pandas.read_csv(file_path)
    # Extract data
    X_encodeur = df.X_Encodeur
    Y_encodeur = df.Y_Encodeur
    X_voulu = df.X_Voulu
    Y_voulu = df.Y_Voulu

    trajectoire_voulue = np.array([X_voulu, Y_voulu])
    trajectoire_obtenue = np.array([X_encodeur, Y_encodeur])

    # Calcul de la distance point par point
    distance = compute_distance(trajectoire_voulue, trajectoire_obtenue)  
    colors = distance
    cmap = plt.get_cmap("viridis")
    # Plot results
    fig, ax = plt.subplots(figsize = (10,8), layout = "constrained")
# 1. Graduations Majeures (tous les 10 cm = 100 mm)
    intervalle_majeur = 100 
    ax.xaxis.set_major_locator(ticker.MultipleLocator(intervalle_majeur))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(intervalle_majeur))
    
    # 2. Graduations Mineures (ex: tous les 1 cm = 10 mm)
    # Vous pouvez changer '10' par '50' (5cm) ou '20' (2cm) selon votre besoin.
    intervalle_mineur = 10 
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(intervalle_mineur))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(intervalle_mineur))
    
    # 3. Style de la grille Majeure (10 cm) : Foncé, continu
    ax.grid(which='major', color='#333333', linestyle='-', linewidth=1.0, alpha=0.8)
    
    # 4. Style de la grille Mineure (1 cm) : Clair, pointillé
    ax.grid(which='minor', color='#CCCCCC', linestyle=':', linewidth=0.5, alpha=0.5)
    
    # S'assurer que les grilles sont dessinées sous les points de données
    ax.set_axisbelow(True)
    # ------------------------------------------
    ax.scatter(X_voulu*1000, Y_voulu*1000, color = 'red', label = 'Commande envoyée')
    sc = ax.scatter(X_encodeur*1000, Y_encodeur*1000, c = distance*1000, cmap = cmap, label = 'Position réelle')
    cbar = plt.colorbar(sc)
    cbar.set_label("Distance d\'erreur [mm]")
    ax.set_title(f"Erreur de position moyenne : {np.mean(distance)*1000:.2f} [mm] et écart-type: {np.std(distance)*1000:.2f} [mm]")
    ax.set_xlabel("Position en X [mm]")
    ax.set_ylabel("Position en Y [mm]")
    ax.legend()
    ax.legend(loc='best')
    ax.set_aspect('equal', adjustable='box')
    # plt.tight_layout() 
    
    # Sauvegarde haute qualité
    fig.savefig(f"Distance_{name}.png", dpi=300, bbox_inches='tight')
    # Optionnel : sauvegarde une version vectorielle pour tes rapports
    fig.savefig(f"Distance_{name}.svg", bbox_inches='tight')
    
    plt.show()

show_trajectories_with_distance(file_path_cercle)
show_trajectories_with_distance(file_path_carre)

aleatoire = extract_data_from_csv(file_path_aleatoire)
show_trajectories(aleatoire)


