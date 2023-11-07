/*
    ================================================================================
    Copyright Centre National de la Recherche Scientifique (CNRS)
        Contributors and copyright holders :

        Xavier Martinez, 2017-2021
        Marc Baaden, 2010-2021
        baaden@smplinux.de
        http://www.baaden.ibpc.fr

        This software is a computer program based on the Unity3D game engine.
        It is part of UnityMol, a general framework whose purpose is to provide
        a prototype for developing molecular graphics and scientific
        visualisation applications. More details about UnityMol are provided at
        the following URL: "http://unitymol.sourceforge.net". Parts of this
        source code are heavily inspired from the advice provided on the Unity3D
        forums and the Internet.

        This program is free software: you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU General Public License for more details.

        You should have received a copy of the GNU General Public License
        along with this program. If not, see <https://www.gnu.org/licenses/>.

        References : 
        If you use this code, please cite the following reference :         
        Z. Lv, A. Tek, F. Da Silva, C. Empereur-mot, M. Chavent and M. Baaden:
        "Game on, Science - how video game technology may help biologists tackle
        visualization challenges" (2013), PLoS ONE 8(3):e57990.
        doi:10.1371/journal.pone.0057990
       
        If you use the HyperBalls visualization metaphor, please also cite the
        following reference : M. Chavent, A. Vanel, A. Tek, B. Levy, S. Robert,
        B. Raffin and M. Baaden: "GPU-accelerated atom and dynamic bond visualization
        using HyperBalls, a unified algorithm for balls, sticks and hyperboloids",
        J. Comput. Chem., 2011, 32, 2924

    Please contact unitymol@gmail.com
    ================================================================================
*/


using UnityEngine;
using System.Collections.Generic;
using System.Linq;

namespace UMol {

public class UnityMolCartoonManager : UnityMolGenericRepresentationManager {

    private List<GameObject> meshesGO;

    public CartoonRepresentation atomRep;

    private bool needUpdate = false;
    // private Material highlightMat;


    /// <summary>
    /// Disables the renderers for all objects managed by the instance of the manager.
    /// </summary>
    public override void DisableRenderers() {
        isEnabled = false;
        if (meshesGO == null)
            return;
        for (int i = 0; i < meshesGO.Count; i++) {
            meshesGO[i].GetComponent<Renderer>().enabled = false;
        }
        // UnityMolMain.getRepresentationManager().UpdateActiveColliders();
    }

    /// <summary>
    /// Enables the renderers for all objects managed by the instance of the manager.
    /// </summary>
    public override void EnableRenderers() {
        isEnabled = true;
        if (meshesGO == null)
            return;
        if (needUpdate) {
            atomRep.recompute();
            needUpdate = false;
        }
        for (int i = 0; i < meshesGO.Count; i++) {
            meshesGO[i].GetComponent<Renderer>().enabled = true;
        }
        // UnityMolMain.getRepresentationManager().UpdateActiveColliders();
    }

    /// <summary>
    /// Initializes this instance of the manager.
    /// </summary>
    public override void Init(SubRepresentation umolRep) {



        if (isInit) {
            return;
        }

        atomRep = (CartoonRepresentation) umolRep.atomRep;
        meshesGO = atomRep.meshesGO;
        // highlightMat = (Material) Resources.Load("Materials/HighlightMaterial");


        isInit = true;
        isEnabled = true;
        areSideChainsOn = true;
        areHydrogensOn = true;
        isBackboneOn = true;
    }

    public override void Clean() {
        GameObject parent = null;
        if (meshesGO.Count != 0) {
            parent = meshesGO[0].transform.parent.gameObject;
        }

        for (int i = 0; i < meshesGO.Count; i++) {
            GameObject.DestroyImmediate(meshesGO[i]);
        }
        if (parent != null) {
            GameObject.DestroyImmediate(parent);
        }

        meshesGO.Clear();
        meshesGO = null;
        atomRep = null;
        isEnabled = false;
        isInit = false;
        // UnityMolMain.getRepresentationManager().UpdateActiveColliders();
    }

    public override void ShowShadows(bool show) {
        if (meshesGO == null)
            return;
        for (int i = 0; i < meshesGO.Count; i++) {
            if (show)
                meshesGO[i].GetComponent<Renderer>().shadowCastingMode = UnityEngine.Rendering.ShadowCastingMode.On;
            else
                meshesGO[i].GetComponent<Renderer>().shadowCastingMode = UnityEngine.Rendering.ShadowCastingMode.Off;
        }
    }

    public override void ShowHydrogens(bool show) {
        areHydrogensOn = show;
    }//Does not really make sense for cartoon

    public override void ShowSideChains(bool show) {
        areSideChainsOn = show;
    }//Does not really make sense for cartoon

    public override void ShowBackbone(bool show) {
        isBackboneOn = show;
    }//Does not really make sense for cartoon

    public void ApplyColors() {
        foreach (GameObject go in meshesGO) {
            go.GetComponent<MeshFilter>().sharedMesh.SetColors(atomRep.meshColors[go]);
        }
    }



    public override void SetColor(Color col, UnityMolSelection sele) {
        List<UnityMolResidue> residues = new List<UnityMolResidue>();

        foreach (UnityMolAtom a in sele.atoms) {
            residues.Add(a.residue);
        }
        var residuesIE = residues.Distinct();
        foreach (UnityMolResidue r in residuesIE) {
            SetColor(col, r, false);
        }
        ApplyColors();
    }

    public void SetColor(Color col, UnityMolResidue r, bool applyNow) {
        List<int> listVertId;
        GameObject curGo = null;

        if (atomRep.residueToGo.TryGetValue(r, out curGo) && atomRep.residueToVert.TryGetValue(r, out listVertId)) {

            List<Color32> colors = atomRep.meshColors[curGo];

            foreach (int c in listVertId) {
                if (c >= 0 && c < colors.Count) {
                    colors[c] = col;
                }
            }

            if (applyNow) {
                ApplyColors();
            }
        }
    }

    public override void SetColor(Color col, UnityMolAtom atom) {//Does not really make sense for atoms
        Debug.LogWarning("Cannot set the color of one atom with the cartoon representation");
    }

    public override void SetColors(Color col, List<UnityMolAtom> atoms) {
        SetColor(col, new UnityMolSelection(atoms, newBonds: null, "tmp"));
    }

    public override void SetColors(List<Color> cols, List<UnityMolAtom> atoms) {
        if (cols.Count == atoms.Count) {
            HashSet<UnityMolResidue> doneRes = new HashSet<UnityMolResidue>();
            for (int i = 0; i < cols.Count; i++) {
                if (!doneRes.Contains(atoms[i].residue)) {
                    SetColor(cols[i], atoms[i].residue, false);
                    doneRes.Add(atoms[i].residue);
                }
            }
            ApplyColors();
        }
        else {
            Debug.LogError("Number of colors should be equal to the number of atoms");
            return;
        }
    }


    public override void SetDepthCueingStart(float v) {
        if (meshesGO == null)
            return;
        foreach (GameObject meshGO in meshesGO) {

            Material[] mats = meshGO.GetComponent<Renderer>().sharedMaterials;
            mats[0].SetFloat("_FogStart", v);
        }
    }

    public override void SetDepthCueingDensity(float v) {
        if (meshesGO == null)
            return;
        foreach (GameObject meshGO in meshesGO) {

            Material[] mats = meshGO.GetComponent<Renderer>().sharedMaterials;
            mats[0].SetFloat("_FogDensity", v);
        }
    }

    public override void EnableDepthCueing() {
        if (meshesGO == null)
            return;
        foreach (GameObject meshGO in meshesGO) {

            Material[] mats = meshGO.GetComponent<Renderer>().sharedMaterials;
            mats[0].SetFloat("_UseFog", 1.0f);
        }
    }

    public override void DisableDepthCueing() {
        if (meshesGO == null)
            return;
        foreach (GameObject meshGO in meshesGO) {

            Material[] mats = meshGO.GetComponent<Renderer>().sharedMaterials;
            mats[0].SetFloat("_UseFog", 0.0f);
        }
    }

    public override void updateWithTrajectory() {

        needUpdate = true;
        bool wasEnabled = true;
        if (meshesGO != null && meshesGO.Count >= 1)
            wasEnabled = meshesGO[0].GetComponent<Renderer>().enabled;

        if (wasEnabled) {
            atomRep.recompute();
            needUpdate = false;
        }
    }

    public override void updateWithModel() {
        needUpdate = true;
        bool wasEnabled = true;
        if (meshesGO != null && meshesGO.Count >= 1)
            wasEnabled = meshesGO[0].GetComponent<Renderer>().enabled;

        if (wasEnabled) {
            atomRep.recompute(isNewModel: true);
            needUpdate = false;
        }
    }

    public override void ShowAtom(UnityMolAtom atom, bool show) {
        Debug.LogWarning("Cannot show/hide one atom with the cartoon representation");
    }

    public override void SetSize(UnityMolAtom atom, float size) {
        Debug.LogWarning("Cannot set the size of one atom with the cartoon representation");
    }

    public override void SetSizes(List<UnityMolAtom> atoms, List<float> sizes) {
        Debug.LogWarning("Cannot set the size of atoms with the cartoon representation");
    }

    public override void SetSizes(List<UnityMolAtom> atoms, float size) {
        Debug.LogWarning("Cannot set the size of atoms with the cartoon representation");
    }
    public override void ResetSize(UnityMolAtom atom) {
        Debug.LogWarning("Cannot set the size of one atom with the cartoon representation");
    }
    public override void ResetSizes() {
        Debug.LogWarning("Cannot set the size of atoms with the cartoon representation");
    }

    public override void ResetColor(UnityMolAtom atom) {
        SetColor(atom.color, atom);
    }
    public override void ResetColors() {

        if (atomRep.savedColors.Count == meshesGO.Count) {
            int i = 0;
            foreach (GameObject go in meshesGO) {
                atomRep.meshColors[go] = atomRep.savedColors[i].ToList();
                i++;
            }
            ApplyColors();
        }
        atomRep.colorationType = colorType.defaultCartoon;
    }

    public override void HighlightRepresentation() {
        // foreach (GameObject meshGO in meshesGO) {

        //     Material[] mats = meshGO.GetComponent<Renderer>().sharedMaterials;

        //     if (mats.Length != 2) {
        //         Material[] newMats = new Material[2];
        //         newMats[0] = mats[0];
        //         newMats[1] = highlightMat;
        //         mats = newMats;
        //         meshGO.GetComponent<Renderer>().sharedMaterials = newMats;
        //     }
        // }

    }


    public override void DeHighlightRepresentation() {

        // foreach (GameObject meshGO in meshesGO) {

        //     Material[] mats = meshGO.GetComponent<Renderer>().sharedMaterials;

        //     if (mats.Length != 1) {
        //         Material[] newMats = new Material[1];
        //         newMats[0] = mats[0];
        //         meshGO.GetComponent<Renderer>().sharedMaterials = newMats;
        //     }
        // }
    }
    public override void SetSmoothness(float val) {
        foreach (GameObject meshGO in meshesGO) {

            Material[] mats = meshGO.GetComponent<Renderer>().sharedMaterials;
            mats[0].SetFloat("_Glossiness", val);
        }
    }
    public override void SetMetal(float val) {
        foreach (GameObject meshGO in meshesGO) {

            Material[] mats = meshGO.GetComponent<Renderer>().sharedMaterials;
            mats[0].SetFloat("_Metallic", val);
        }
    }
    public override UnityMolRepresentationParameters Save() {
        UnityMolRepresentationParameters res = new UnityMolRepresentationParameters();
        res.repT.atomType = AtomType.cartoon;
        res.repT.bondType = BondType.nobond;
        res.colorationType = atomRep.colorationType;

        if (res.colorationType == colorType.custom) {
            res.colorPerAtom = new Dictionary<UnityMolAtom, Color32>(atomRep.savedAtoms.Count);
            foreach (UnityMolAtom a in atomRep.savedAtoms) {
                List<int> listVertId;
                GameObject curGo = null;
                if (atomRep.residueToGo.TryGetValue(a.residue, out curGo) && atomRep.residueToVert.TryGetValue(a.residue, out listVertId)) {
                    List<Color32> colors = atomRep.meshColors[curGo];
                    res.colorPerAtom[a] = colors[listVertId[0]];
                }
            }
        }
        else if (res.colorationType == colorType.full) { //Get color of first atom/residue
            foreach (UnityMolAtom a in atomRep.savedAtoms) {
                List<int> listVertId;
                GameObject curGo = null;
                if (atomRep.residueToGo.TryGetValue(a.residue, out curGo) && atomRep.residueToVert.TryGetValue(a.residue, out listVertId)) {
                    List<Color32> colors = atomRep.meshColors[curGo];
                    res.fullColor = colors[0];
                    break;
                }
            }
        }
        else if (res.colorationType == colorType.bfactor) {
            res.bfactorStartColor = atomRep.bfactorStartCol;
            res.bfactorEndColor = atomRep.bfactorEndCol;
        }
        if (meshesGO != null && meshesGO.Count >= 1) {
            res.smoothness = meshesGO[0].GetComponent<Renderer>().sharedMaterials[0].GetFloat("_Glossiness");
            res.metal = meshesGO[0].GetComponent<Renderer>().sharedMaterials[0].GetFloat("_Metallic");
            res.shadow = (meshesGO[0].GetComponent<Renderer>().shadowCastingMode == UnityEngine.Rendering.ShadowCastingMode.On);
        }

        return res;
    }

    public override void Restore(UnityMolRepresentationParameters savedParams) {

        if (savedParams.repT.atomType == AtomType.cartoon && savedParams.repT.bondType == BondType.nobond) {
            if (savedParams.colorationType == colorType.full) {
                SetColor(savedParams.fullColor, atomRep.selection);
                atomRep.colorationType = colorType.full;
            }
            else if (savedParams.colorationType == colorType.custom) {
                foreach (UnityMolAtom a in atomRep.selection.atoms) {
                    if (savedParams.colorPerAtom.ContainsKey(a)) {
                        SetColor(savedParams.colorPerAtom[a], a.residue, false);
                    }
                }
                ApplyColors();
            }
            else if (savedParams.colorationType == colorType.defaultCartoon) {
                //Do nothing !
            }
            else if (savedParams.colorationType == colorType.res) {
                colorByRes(atomRep.selection);
            }
            else if (savedParams.colorationType == colorType.chain) {
                colorByChain(atomRep.selection);
            }
            else if (savedParams.colorationType == colorType.hydro) {
                colorByHydro(atomRep.selection);
            }
            else if (savedParams.colorationType == colorType.seq) {
                colorBySequence(atomRep.selection);
            }
            else if (savedParams.colorationType == colorType.charge) {
                colorByCharge(atomRep.selection);
            }
            else if (savedParams.colorationType == colorType.restype) {
                colorByResType(atomRep.selection);
            }
            else if (savedParams.colorationType == colorType.rescharge) {
                colorByResCharge(atomRep.selection);
            }
            else if (savedParams.colorationType == colorType.bfactor) {
                colorByBfactor(atomRep.selection, savedParams.bfactorStartColor, savedParams.bfactorEndColor);
            }

            SetMetal(savedParams.metal);
            SetSmoothness(savedParams.smoothness);
            ShowShadows(savedParams.shadow);
            atomRep.colorationType = savedParams.colorationType;
        }
        else {
            Debug.LogError("Could not restore representation parameteres");
        }
    }
}
}