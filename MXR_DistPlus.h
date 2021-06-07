/*
 * MXR_DistPlus.h
 *
 *  Created on: May 8, 2021
 *      Author: Sam S. AME
 */

#ifndef MXR_DISTPLUS_H_
#define MXR_DISTPLUS_H_

#include "wdf_utils.h"

class mxrDistPlus1 : public wdfTree
{
private:
	ResistiveVoltageSource* Vres;
	ResistiveVoltageSource* Vb; //encompasses R2

	Resistor* R1;
	Resistor* R3;
	Resistor* R4;

	Resistor* ResDist; //distortion potentiometer
	float ResDistPrev;

	Capacitor* C1;
	Capacitor* C2;
	Capacitor* C3;
	Capacitor* C4;

	Series* S1;
	Series* S2;
	Series* S3;
	Series* S4;

	Parallel* P1;
	Parallel* P2;

	RootRtypeAdaptor* R;
	wdfNode** R_subTreeNodes;

public:

	mxrDistPlus1()
	{
		ResDistPrev = 0.0;

		//Port A
		Vres = new ResistiveVoltageSource(1.0);
		Vb = new ResistiveVoltageSource(1e6);
		Vb->setVoltage(4.5);
		R1 = new Resistor(10000.0);
		C1 = new Capacitor(1.0e-9,sr,alpha);
		C2 = new Capacitor(10.0e-9,sr,alpha);
		S3 = new Series(C2,R1);
		P1 = new Parallel(Vres,C1);
		S4 = new Series(S3,P1);
		P2 = new Parallel(S4,Vb);

		//Port B
		ResDist = new Resistor(1e6);
		R3 = new Resistor(4700.0);
		C3 = new Capacitor(47.0e-9,sr,alpha);
		S1 = new Series(ResDist,R3);
		S2 = new Series(S1,C3);

		//Port C
		R4 = new Resistor(1e6);

		//Port D
		C4 = new Capacitor(1.0e-6,sr,alpha);

		//Add ports to subtree nodes
		R_subTreeNodes = new wdfNode*[4];
		R_subTreeNodes[0] = P2;
		R_subTreeNodes[1] = S2;
		R_subTreeNodes[2] = R4;
		R_subTreeNodes[3] = C4;

		R = new RootRtypeAdaptor(4,R_subTreeNodes);

		setSMatrixData();
	}

	~mxrDistPlus1()
	{
		delete R1; delete R3; delete R4; delete ResDist;
		delete C1; delete C2; delete C3; delete C4;
		delete S1; delete S2; delete S3; delete S4;
		delete P1; delete P2;
		delete Vres; delete Vb;
		delete R;
		delete R_subTreeNodes;
	}

	void setSMatrixData()
	{
		float Ra = R_subTreeNodes[0]->Rp;
		float Rb = R_subTreeNodes[1]->Rp;
		float Rc = R_subTreeNodes[2]->Rp;
		float Rd = R_subTreeNodes[3]->Rp;

		double arr[4][4] =
		{
				{                                     1,                                0,                                0,  0},
				{                (200*Rb)/(101*Rb + Rc),       1 - (202*Rb)/(101*Rb + Rc),             (2*Rb)/(101*Rb + Rc),  0},
				{               -(200*Rc)/(101*Rb + Rc),           (202*Rc)/(101*Rb + Rc),         1 - (2*Rc)/(101*Rb + Rc),  0},
				{(200*Rd*(Rb + Rc))/(101*Rb*Rd + Rc*Rd), -(200*Rc*Rd)/(101*Rb*Rd + Rc*Rd), -(200*Rb*Rd)/(101*Rb*Rd + Rc*Rd), -1}
		};

		for(int i = 0; i < 4; i++)
		{
			for(int j = 0; j < 4; j++)
			{
				R->getSMatrixData()[i][j] = static_cast<float>(arr[i][j]);
			}
		}
	}

	float processSample(float inSamp)
	{
		Vres->setVoltage(inSamp);

		R->calcIncidentWave(0.0);
		R->calcReflectedWave();

		return C4->getPortVoltage();
	}

	void setParams(float rDist)
	{
		if(rDist != ResDistPrev)
		{
			ResDist->setResistance(rDist);
			ResDistPrev = rDist;
			setSMatrixData();
		}
	}

};

class mxrDistPlus2 : public wdfTree
{
private:
	ResistiveVoltageSource* Vres; //encompasses R5 = 10k
	Resistor* ResOp; //output potentiometer
	Resistor* ResOm;
	DiodePair* DP;
	Capacitor* C5;
	Series* S5;
	Parallel* P3;
	Parallel* P4;

public:

	mxrDistPlus2()
	{
		Vres = new ResistiveVoltageSource(10000.0);

		ResOp = new Resistor(5000.0);
		ResOm = new Resistor(5000.0);

		C5 = new Capacitor(1.0e-9,sr,alpha);

		DP = new DiodePair(d_shockley);

		S5 = new Series(ResOp, ResOm);

		P3 = new Parallel(S5,C5);
		P4 = new Parallel(P3,Vres);

		P4->connectToNode(DP);
	}

	void setSMatrixData(){}

	float processSample(float inSamp)
	{
		Vres->setVoltage(inSamp);

		DP->calcIncidentWave(P4->calcReflectedWave());
		P4->calcIncidentWave(DP->calcReflectedWave());

		return ResOm->getPortVoltage();
	}

	void setParams(float alpha)
	{
		ResOp->setResistance(alpha*10000.0);
		ResOm->setResistance((1.0-alpha)*10000.0);
	}
};

#endif /* MXR_DISTPLUS_H_ */
