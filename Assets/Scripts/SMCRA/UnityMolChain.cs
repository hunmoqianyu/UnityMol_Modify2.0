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
/// <summary>
/// Part of the SMCRA data structure, UnityMolChain stores the residues of the structure as a list of UnityMolResidue
/// A reference to the model it belongs is provided
/// </summary>
public class UnityMolChain {

	/// <summary>
	/// Store all the residues of the chain based on their ids
	/// </summary>
	public Dictionary<int, UnityMolResidue> residues;

	/// <summary>
	/// Reference to the model the chain belongs to
	/// </summary>
	public UnityMolModel model;

	/// <summary>
	/// Name of the chain
	/// </summary>
	public string name;


	/// <summary>
	/// UnityMolChain constructor taking a list of residues as arg
	/// </summary>
	public UnityMolChain(List<UnityMolResidue> _residues, string _name) {
		residues = new Dictionary<int, UnityMolResidue>();
		AddResidues(_residues);
		name = _name;
	}

	/// <summary>
	/// UnityMolChain constructor taking a residue as arg
	/// </summary>
	public UnityMolChain(UnityMolResidue _residue, string _name) {
		residues = new Dictionary<int, UnityMolResidue>();
		residues[_residue.id] = _residue;
		name = _name;
	}

	/// <summary>
	/// Add a list of residues to the stored residues
	/// </summary>
	public void AddResidues(List<UnityMolResidue> newResidues) {
		foreach (UnityMolResidue r in newResidues) {
			residues[r.id] = r;
		}
	}

	/// <summary>
	/// Add a dictionary of residues to the stored residues
	/// </summary>
	public void AddResidues(Dictionary<int, UnityMolResidue> newResidues) {
		foreach (UnityMolResidue r in newResidues.Values) {
			residues[r.id] = r;
		}
	}

	public List<UnityMolAtom> allAtoms {
		get { return ToAtomList(); }
	}

	public List<UnityMolAtom> ToAtomList() {
		List<UnityMolAtom> res = new List<UnityMolAtom>();

		foreach (UnityMolResidue r in residues.Values) {
			// res.AddRange(r.allAtoms);
			foreach (UnityMolAtom a in r.atoms.Values) {
				res.Add(a);
			}
		}
		return res;
	}

	public UnityMolSelection ToSelection(bool doBonds = true) {
		List<UnityMolAtom> selectedAtoms = ToAtomList();
		string selectionMDA = model.structure.uniqueName + " and chain " + name;

		if (doBonds) {
			return new UnityMolSelection(selectedAtoms, name, selectionMDA);
		}
		return new UnityMolSelection(selectedAtoms, newBonds: null, name, selectionMDA);
	}
	/**
	* 氨基酸之间脱水缩合，res_N是待脱水的氨基所在位置，res_C是待脱水的羧基所在位置
	* 其中羧基COOH所在的氨基酸必须是独立的，没有和其他氨基酸形成链
	* 对于氨基NH2没有要求
	*/
	public void DehydrationAndCondensation(int res_N,int res_C){
		Vector3 H1_Position=new Vector3(0,0,0);
		Vector3 C_Position=new Vector3(0,0,0);
		Vector3 H2_Position=new Vector3(0,0,0);
		Vector3 N_Position=new Vector3(0,0,0);
		Vector3 O_Position=new Vector3(0,0,0);
		Vector3 CA_Position=new Vector3(0,0,0);

		// 记录C,N原子，用于生成化学键 
		UnityMolAtom C=null;
		UnityMolAtom N=null;
		UnityMolAtom H1=null;
		UnityMolAtom OXT=null;
		UnityMolAtom HXT=null;
		// 移除原子，记录对应原子坐标信息
		foreach(UnityMolAtom r in residues[res_N].allAtoms){
			if(r.name=="H1"){
				H1=r;
				H1_Position=r.oriPosition;
				residues[res_N].atoms.Remove(r.name);//移除H原子
			}
			else if(r.name=="C"){
				CA_Position=r.oriPosition;
			}
			else if(r.name=="H2"){
				H2_Position=r.oriPosition;
			}
			else if(r.name=="N"){
				N=r;
				N_Position=r.oriPosition;
			}
		}
		foreach(UnityMolAtom r in residues[res_C].allAtoms){
			if(r.name=="C"){
				C=r;
				C_Position=r.oriPosition;
			}
			else if(r.name=="O"){
				O_Position=r.oriPosition;
			}
			else if(r.name=="OXT"){
				OXT=r;
				residues[1].atoms.Remove(r.name);//移除O,H原子
			}
			else if(r.name=="HXT"){
				HXT=r;
				residues[1].atoms.Remove(r.name);
			}
		}
		//增加C,N的化学键
		model.bonds.Add(C,N);
		//移除无用的化学键
		model.bonds.Remove(N,H1);
		model.bonds.Remove(C,OXT);
		model.bonds.Remove(OXT,HXT);




		// 调整原子坐标，使得肽键在同一个平面上
		H1_Position=(H1_Position-N_Position)*1.35f+N_Position;//将键长扩展成1.35倍
		Vector3 offset=H1_Position-C_Position;//计算平移偏移
		C_Position+=offset;//更新位置信息，便于计算旋转角和旋转轴
		O_Position+=offset;
		Vector3 normal=Vector3.Cross(H2_Position-N_Position,C_Position-N_Position).normalized;//计算法向量
		float rotationAngle=90f-Vector3.Angle(normal,O_Position-C_Position);//计算旋转角
		Vector3 rotationAxis=Vector3.Cross(normal,O_Position-C_Position);//计算旋转轴
		Quaternion rotation =Quaternion.AngleAxis(rotationAngle,rotationAxis);//计算旋转变换
		foreach(UnityMolAtom r in residues[res_C].allAtoms){
			r.oriPosition=r.oriPosition+offset;//平移
			r.oriPosition=rotation*(r.oriPosition-C_Position)+C_Position;//以C为中心进行旋转
		}
		//优化，继续旋转，以C为中心，旋转轴为计算出来的法向量normal,
        float best=0;//最好的旋转角度
		float max_distance=0;
		float length=(N_Position-CA_Position).magnitude;//C,N向量的模长
		for(float angle=180f*(residues.Count%2);angle<180f*(residues.Count%2+1); angle+=30f){//交替旋转
			float distance=0;
			Quaternion trial_r=Quaternion.AngleAxis(angle,normal);
			foreach(UnityMolAtom r in residues[res_C].allAtoms){
			    Vector3 temp=trial_r*(r.oriPosition-C_Position)+C_Position-CA_Position;//计算原子与N原子的向量
				distance+=(Vector3.Dot(temp,N_Position-CA_Position)/length);
				
		    }	
			if(distance>max_distance){
				max_distance=distance;
				best=angle;
			}
		
		}
		rotation =Quaternion.AngleAxis(best,normal);//计算旋转变换
		foreach(UnityMolAtom r in residues[res_C].allAtoms){
			r.oriPosition=rotation*(r.oriPosition-C_Position)+C_Position;//以C为中心进行旋转
		}

	}

}
}