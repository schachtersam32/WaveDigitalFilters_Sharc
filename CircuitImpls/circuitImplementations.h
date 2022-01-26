/*
 * circuitImplementation.h
 *
 *  Created on: Apr 6, 2021
 *      Author: Sam S. AME
 */

#ifndef CIRCUITIMPLEMENTATION_H_
#define CIRCUITIMPLEMENTATION_H_

#include "wdf_utils.h"

#define FOUR_PI_SQ 4*PI*PI

class TwoPoleButterworthLowPass : public wdfTree
{
private:
	Resistor* R1;
	Resistor* RL;
	Capacitor* C1;
	Inductor* L1;

	Series* S1;
	Series* S2;
	Parallel* P1;

	IdealVoltageSource* Vin;

	float fc; //cutoff frequency

public:
	TwoPoleButterworthLowPass() : fc(0.0)
	{
		R1 = new Resistor(50.0);
		RL = new Resistor(1000000.0);
		C1 = new Capacitor(47e-6,sr,alpha);
		L1 = new Inductor(110e-6,sr,alpha);

		Vin = new IdealVoltageSource;

		S1 = new Series(RL,L1);
		P1 = new Parallel(S1,C1);
		S2 = new Series(R1,P1);

		S2->connectToNode(Vin);
	}

	~TwoPoleButterworthLowPass()
	{
		delete R1; delete RL;
		delete C1; delete L1;
		delete S1; delete S2; delete P1;
		delete Vin;
	}

	float processSample(float inSamp)
	{
		Vin->setVoltage(inSamp);

		Vin->calcIncidentWave(S2->calcReflectedWave());
		float output = RL->getPortVoltage();
		S2->calcIncidentWave(Vin->calcReflectedWave());

		return output;
	}

	void setParams(float cutoffFreq)
	{
		if(cutoffFreq != fc)
		{
			fc = cutoffFreq;
			L1->setInductance(R1->getResistance()*RL->getResistance()*C1->getCapacitance()*fc*fc*FOUR_PI_SQ);
		}
	}

};

class PassiveToneControl : public wdfTree
{
private:

	IdealVoltageSource* Vin;
	Resistor* Res1;
	Resistor* Res2;
	Resistor* Res3p;
	Resistor* Res3m;
	Capacitor* C1;
	Capacitor* C2;

	Series* S1;

	RtypeAdaptor* R1;
	wdfNode** R_subtreeNodes;

	float pot_val;
	float R3;
public:

	PassiveToneControl()
	{
		pot_val = 0.5; R3 = 250000;

		Vin = new IdealVoltageSource; // port A

		Res1 = new Resistor(47000); // port C
		Res2 = new Resistor(47000); // port E
		Res3p = new Resistor(pot_val*R3);
		Res3m = new Resistor((1.0-pot_val)*R3);
		S1 = new Series(Res3p,Res3m); // port F
		C1 = new Capacitor(4.7e-6,sr,alpha); // port B
		C2 = new Capacitor(0.022e-6,sr,alpha); // port D

		R_subtreeNodes = new wdfNode*[6];
		R_subtreeNodes[0] = Vin;
		R_subtreeNodes[1] = C1;
		R_subtreeNodes[2] = Res1;
		R_subtreeNodes[3] = C2;
		R_subtreeNodes[4] = Res2;
		R_subtreeNodes[5] = S1;

		R1 = new RtypeAdaptor(6,R_subtreeNodes);
	}

	void setSMatrixData()
	{
		float arr[6][6] =
		{
			 {0,    0.5454,    0.4600,    0.5400,    0.4546,   -0.0855},
			 {0.0001,    0.9999,   -0.0001,    0.0000 ,   0.0000,    0.0000},
			 {0.9999,   -1.4545,   -0.5400,    0.5400 ,   0.4545,   -0.0855},
			 {0.0118,    0.0027,    0.0054,    0.9827 ,  -0.0145,    0.0027},
			 {0.9882,    0.5427,    0.4545,   -1.4427 ,  -0.5309,   -0.0882},
			 {-0.9881,    1.4572,   -0.4546,    1.4427,   -0.4691,   -0.9118}
		};

		for(int i = 0; i < 4; i++)
		{
			for(int j = 0; j < 4; j++)
			{
				R1->getSMatrixData()[i][j] = arr[i][j];
			}
		}

	}

	float processSample(float inSamp)
	{
		Vin->setVoltage(inSamp);

		R1->calcIncidentWave(0.0);
		R1->calcReflectedWave();

		return Res3m->getPortVoltage();
	}

	void setParams(float val)
	{
		if(val != pot_val)
		{
			pot_val = val;
			Res3p->setResistance(pot_val*R3);
			Res3m->setResistance((1-pot_val)*R3);
		}
	}

};

class marshallDiodeClipper : public wdfTree
{
private:
	DiodePair* DP1;
	Capacitor* C3;
	ResistiveCurrentSource* Is1;
	Parallel* P3;

public:
	marshallDiodeClipper()
	{
		setAlpha(0.0);

		C3 = new Capacitor(47.0e-12, sr, alpha);
		Is1 = new ResistiveCurrentSource(220000.0);
		DP1 = new DiodePair(d_NSCW100);
		P3 = new Parallel(Is1,C3);

		P3->connectToNode(DP1);
	}

	~marshallDiodeClipper()
	{
		delete C3;
		delete Is1;
		delete DP1;
		delete P3;
	}

	float processSample(float inCurrent)
	{
		Is1->setCurrent(inCurrent);

		DP1->calcIncidentWave(P3->calcReflectedWave());
		P3->calcIncidentWave(DP1->calcReflectedWave());

		return Is1->getPortVoltage();
	}

	void setSMatrixData() {}

	void setParams(float R_val)
	{
		Is1->setResistance(R_val);
	}

};

class marshallBuffer : public wdfTree
{
private:
	Resistor* R1;
	Resistor* R2;
	Resistor* R_pot;
	Capacitor* C1;
	Capacitor* C2;
	IdealVoltageSource* Vin;
	Series* S1;
	Series* S2;
	Parallel* P1;
	Parallel* P2;
public:
	marshallBuffer()
	{
		setAlpha(0.0);

		R1 = new Resistor(22000.0);
		R2 = new Resistor(12000.0);
		R_pot = new Resistor(1.0);

		C1 = new Capacitor(47e-9,sr,alpha);
		C2 = new Capacitor(1e-9,sr,alpha);

		Vin = new IdealVoltageSource;

		P2 = new Parallel(R_pot,C2);
		S2 = new Series(P2,R2);
		P1 = new Parallel(S2,R1);
		S1 = new Series(P1,C1);

		S1->connectToNode(Vin);
	}

	~marshallBuffer()
	{
		delete R1;
		delete R2;
		delete R_pot;
		delete C1;
		delete C2;
		delete Vin;
		delete P2;
		delete P1;
		delete S2;
		delete S1;
	}

	float processSample(float inSamp)
	{
		Vin->setVoltage(inSamp);

		Vin->calcIncidentWave(S1->calcReflectedWave());
		S1->calcIncidentWave(Vin->calcReflectedWave());

		return R2->getPortCurrent();
	}

	void setSMatrixData() {}

	void setParams(float R_val)
	{
		R_pot->setResistance(R_val);
	}

};

class switchableAttenuator : wdfTree
{
private:
	Switch* SW1;
	Resistor* R1;
	Resistor* R2;
	Parallel* P1;
	Series* S1;
	ResistiveVoltageSource* Vs;
public:
	switchableAttenuator()
	{

		R1 = new Resistor(250000.0);
		R2 = new Resistor(250000.0);
		Vs = new ResistiveVoltageSource(1.0);
		SW1 = new Switch;

		S1 = new Series(Vs,R2);
		P1 = new Parallel(S1,R1);

		P1->connectToNode(Vs);
	}

	~switchableAttenuator(){}

	float processSample(float inSamp)
	{
		Vs->setVoltage(inSamp);

		SW1->calcIncidentWave(P1->calcReflectedWave());
		P1->calcIncidentWave(SW1->calcReflectedWave());

		return R2->getPortVoltage();
	}

	void setParams(bool swClose)
	{
		SW1->closeSwitch(swClose);
	}

};

class DiodeClipper : wdfTree
{
private:
	Capacitor* C1;
	Capacitor* C2;
	DiodePair* DP;
	Series* S1;
	Parallel* P1;
	ResistiveVoltageSource* Vs;

public:
	DiodeClipper()
	{
		setAlpha(0.0);

		C1 = new Capacitor(0.01e-6,sr,alpha);
		C2 = new Capacitor(0.47e-6,sr,alpha);
		Vs = new ResistiveVoltageSource(2200.0);
		DP = new DiodePair(d_1N4148);

		S1 = new Series(Vs,C2);
		P1 = new Parallel(S1,C1);

		P1->connectToNode(DP);
	}

	~DiodeClipper()
	{
		delete C1; delete C2;
		delete DP; delete Vs;
		delete S1; delete P1;
	}
	void setParams(float alpha)
	{
		C1->setCapacitance(0.01e-6 * alpha);
	}

	float processSample(float inSamp)
	{
		Vs->setVoltage(inSamp);

		DP->calcIncidentWave(P1->calcReflectedWave());
		float outSamp = DP->getPortVoltage();
		P1->calcIncidentWave(DP->calcReflectedWave());

		return outSamp;
	}
};



#endif /* CIRCUITIMPLEMENTATION_H_ */
