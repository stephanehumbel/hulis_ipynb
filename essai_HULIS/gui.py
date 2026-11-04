# ==================== interface essai_HULIS ====================
import ipywidgets as w
from IPython.display import display, clear_output
import numpy as np
import molecule, huckel, huckel_analysis, pipeline

state = {"mol": None, "res": None, "lewis": []}      # mémoire de l'interface

def todo(nom):
    print(f"⚠️ TODO : implémenter « {nom} »")

# ---------- widgets ----------
name_box  = w.Text(placeholder="nom PubChem (ex. fulvene)", layout=w.Layout(width='180px'))
go_name   = w.Button(description="PubChem", button_style='primary')
load_xyz  = w.Text(placeholder="fichier .xyz", layout=w.Layout(width='180px'))
go_xyz    = w.Button(description="Charger XYZ", button_style='primary')

gen_all   = w.Button(description="generate all")     # toutes les Lewis auto
b_create  = w.Button(description="create")          # créer UNE Lewis à la main
b_results = w.Button(description="results")         # recouvrements ⟨φi|φj⟩
b_erase1  = w.Button(description="erase 1")          # effacer la dernière
b_erase_all= w.Button(description="erase all")

viz_orb   = w.Button(description="Orbitales π")     # dessiner HOMO/LUMO etc.
report_b  = w.Button(description="Rapport texte")

out = w.Output()

# ---------- actions ----------
def _ensure_mol():
    if state["mol"] is None:
        raise ValueError("Charge d'abord une molécule (PubChem ou XYZ)")
    return state["mol"]

def on_go_name(_):
    with out:
        clear_output()
        try:
            state["mol"], state["res"] = pipeline.mol_from_name(name_box.value), None
            print(f"✓ {name_box.value} chargée :", len(state['mol'].atoms), "atomes")
            print(state["mol"].summary())
        except Exception as ex:
            print("✗", ex)

def on_go_xyz(_):
    with out:
        clear_output()
        try:
            state["mol"], state["res"] = pipeline.mol_from_xyz(load_xyz.value), None
            print(f"✓ {load_xyz.value} chargée")
            print(state["mol"].summary())
        except Exception as ex:
            print("✗", ex)

def on_report(_):
    with out:
        clear_output()
        mol = _ensure_mol()
        res = state["res"] or (state.__setitem__("res", huckel.solve_huckel(mol)) or state["res"])
        print(huckel_analysis.format_report(mol, res)[0])

def on_gen_all(_):
    with out: clear_output(); _ensure_mol(); todo("generate all (Lewis automatiques)")

def on_create(_):
    with out: clear_output(); _ensure_mol(); todo("create (Lewis à la main)")

def on_results(_):
    with out: clear_output(); _ensure_mol(); todo("results (⟨φi|φj⟩)")

def on_erase1(_):
    with out: clear_output(); _ensure_mol(); todo("erase 1")

def on_erase_all(_):
    with out: clear_output(); state["lewis"].clear(); print("✓ toutes les Lewis effacées")

def on_viz(_):
    with out: clear_output(); _ensure_mol(); todo("dessin des orbitales (voir hulis.ipynb)")

go_name.on_click(on_go_name);   go_xyz.on_click(on_go_xyz)
gen_all.on_click(on_gen_all);   b_create.on_click(on_create)
b_results.on_click(on_results); b_erase1.on_click(on_erase1)
b_erase_all.on_click(on_erase_all)
viz_orb.on_click(on_viz);       report_b.on_click(on_report)

# ---------- mise en page ----------
col_src  = w.VBox([w.Label("Source"), name_box, go_name, load_xyz, go_xyz])
col_lewis= w.VBox([w.Label("Structures de Lewis"), gen_all, b_create, b_results, b_erase1, b_erase_all])
col_viz  = w.VBox([w.Label("Visualisation"), viz_orb, report_b])
panel    = w.HBox([col_src, col_lewis, col_viz])
display(panel, out)

