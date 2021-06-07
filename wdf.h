/*
 * wdf.h
 *
 *  Created on: Apr 6, 2021
 *      Author: Sam S. AME
 */

#ifndef WDF_H_
#define WDF_H_

#include <matrix.h>
#include <vector>
#include "wdf_utils.h"
#include <stddef.h>
#include <stdlib.h>
#include "common/audio_system_config.h"


class wdfTree;
class wdfNode;
class wdfPort;

class wdfPort
{
public:
	float a; //incident wave
	float b; //reflected wave
	float Rp; //port impedance

public:
	wdfPort(std::string portType) : portType(portType), a(0.0), b(0.0), Rp(0.0) {}
	virtual ~wdfPort() {}

	virtual void calcImpedance() = 0; //calculates the port impedance
	virtual void propagateImpedance() = 0; //propagate port impedance up the tree

	virtual void calcIncidentWave(float x) = 0; //find the wave leaving the port
	virtual float calcReflectedWave() = 0; //find the wave entering the port

	float getPortVoltage() //convert from wave to kirchhoff variables
	{
		return (a + b) / 2.0;
	}
	float getPortCurrent()
	{
		return (a - b) / (2.0 * Rp);
	}

private:
	std::string portType;
};

class wdfNode : public wdfPort
{
public:
	wdfNode(std::string type) : wdfPort(type), connectedNode(NULL) {} //initialize a node with a port
	virtual ~wdfNode() {}

	void connectToNode(wdfPort* node) { connectedNode = node; } //connect ports together

	void propagateImpedance()
	{
		calcImpedance();

		if (connectedNode != NULL)
			connectedNode->propagateImpedance();
	}

protected:
	wdfPort* connectedNode;

};

class wdfTree //a special instance of this class is created for each circuit to be run
{
public:
	wdfTree() : alpha(1.0), sr(AUDIO_SAMPLE_RATE) {}
	virtual ~wdfTree() {}

	virtual float processSample(float inSamp) = 0; //processes tree sample by sample
	virtual void setSMatrixData() {}

	void setSampleRate(float SR) { sr = AUDIO_SAMPLE_RATE; } //set system samplerate
	void setAlpha(float a) { alpha = a; } //set circuit alpha parameter
protected:
	float sr;
	float alpha;
};

class Resistor : public wdfNode
{
public:
	Resistor(float R) : wdfNode("Resistor"), R_val(R)
	{
		calcImpedance();
	}

	~Resistor() {}

	void setResistance(float R)
	{
		if (R == R_val)
			return;
		R_val = R;
		propagateImpedance();
	}

	float getResistance() {return R_val;}

	void calcImpedance()
	{
		Rp = R_val;
	}

	void calcIncidentWave(float downWave)
	{
		a = downWave;
	}

	float calcReflectedWave()
	{
		b = 0.0;
		return b;
	}

protected:
	float R_val;
};

class RootResistor : public wdfNode
{
public:
	RootResistor(float R) : wdfNode("RootResistor"), R_val(R)
	{
		calcImpedance();
	}

	~RootResistor() {}

	void setResistance(float R)
	{
		if (R == R_val)
			return;
		R_val = R;
		propagateImpedance();
	}

	float getResistance() {return R_val;}

	void calcImpedance() {}

	void calcIncidentWave(float downWave)
	{
		a = downWave;
	}

	float calcReflectedWave()
	{
		b = (R_val - Rp) / (R_val + Rp);
		return b;
	}

protected:
	float R_val;
};

class Capacitor : public wdfNode
{
public:
	Capacitor(float C, float fs, float alpha) : wdfNode("Capacitor"),
		C_val(C), prev_a(0.0), alpha(alpha), fs(fs), b_coef((1 - alpha) / 2.0), a_coef((1 + alpha) / 2.0)
	{
		calcImpedance();
	}

	~Capacitor() {}

	void setCapacitance(float C)
	{
		if (C == C_val)
			return;

		C_val = C;
		propagateImpedance();
	}

	float getCapacitance() {return C_val;}

	void calcImpedance()
	{
		Rp = 1 / (C_val * fs * (1 + alpha));
	}

	void calcIncidentWave(float downSamp)
	{
		a = downSamp;
		prev_a = a;
	}

	float calcReflectedWave()
	{
		b = b_coef * b + a_coef * prev_a;
		return b;
	}

protected:
	float C_val;
private:
	float prev_a; //storage values for incident and reflected waves
	float fs; // sample rate needed for adaptation
	float alpha; //alpha parameter needed for discretization
	const float b_coef;
	const float a_coef;
};

class RootCapacitor : public wdfNode
{
public:
	RootCapacitor(float C, float fs, float alpha) : wdfNode("RootCapacitor"),
		C_val(C), prev_a(0.0), alpha(alpha), fs(fs)
	{
		calcImpedance();
		calcCoefficients();
	}

	~RootCapacitor() {}

	void setCapacitance(float C)
	{
		if (C == C_val)
			return;

		C_val = C;
		//propagateImpedance();
		calcCoefficients();
	}

	float getCapacitance() {return C_val;}

	void calcCoefficients()
	{
		float T = 1.0 / fs;
		float den = 1.0 / (Rp * C_val * (1 + alpha) + T);
		b_coef = (Rp * C_val * (1 + alpha) - T * alpha) * den;
		a_coef = (-Rp * C_val * (1 + alpha) - T) * den;
		a1_coef = (Rp * C_val * (1 + alpha) + T * alpha) * den;
	}

	void calcImpedance() {}

	void calcIncidentWave(float downSamp)
	{
		a = downSamp;
		prev_a = a;
	}

	float calcReflectedWave()
	{
		b = b_coef * b + a1_coef * prev_a + a_coef * a;
		return b;
	}

protected:
	float C_val;
private:
	float prev_a; //storage values for incident and reflected waves
	float fs; // sample rate needed for adaptation
	float alpha; //alpha parameter needed for discretization
	float b_coef;
	float a1_coef;
	float a_coef;
};

class Inductor : public wdfNode
{
public:
	Inductor(float L, float fs, float alpha) : wdfNode("Inductor"),
		L_val(L), prev_a(0.0), alpha(alpha), fs(fs), b_coef((1 - alpha) / 2.0), a_coef((1 + alpha) / 2.0)
	{
		calcImpedance();
	}

	~Inductor() {}

	void setInductance(float L)
	{
		if (L == L_val)
			return;

		L_val = L;
		propagateImpedance();
	}

	float getInductance() {return L_val;}

	void calcImpedance()
	{
		Rp = (L_val * fs * (1.0 + alpha));
	}

	void calcIncidentWave(float downSamp)
	{
		a = downSamp;
		prev_a = a;
	}

	float calcReflectedWave()
	{
		b = b_coef * b - a_coef * prev_a;
		//		prev_b = b;
		return b;
	}


protected:
	float L_val;
private:
	float prev_a, prev_b; //storage values for incident and reflected waves
	float fs; // sample rate needed for adaptation
	float alpha; //alpha parameter needed for discretization
	const float b_coef;
	const float a_coef;
};

class RootInductor : public wdfNode
{
public:
	RootInductor(float L, float fs, float alpha) : wdfNode("Inductor"),
		L_val(L), prev_a(0.0), alpha(alpha), fs(fs)
	{
		calcImpedance();
		calcCoefficients();
	}

	~RootInductor() {}

	void setInductance(float L)
	{
		if (L == L_val)
			return;

		L_val = L;
		propagateImpedance();
	}

	float getInductance() {return L_val;}

	void calcCoefficients()
	{
		float T = 1.0 / fs;
		float den = 1.0 / (Rp * T + L_val * (1 + alpha));
		b_coef = L_val * (1 + alpha) - Rp * T * alpha;
		a1_coef = L_val * (1 + alpha) - Rp * T;
		a_coef = L_val * (1 + alpha) + Rp * T * alpha;
	}

	void calcImpedance()
	{
		Rp = (L_val * fs * (1.0 + alpha));
	}

	void calcIncidentWave(float downSamp)
	{
		a = downSamp;
		prev_a = a;
	}

	float calcReflectedWave()
	{
		b = b_coef * b + a1_coef * prev_a + a_coef * a;
		return b;
	}


protected:
	float L_val;
private:
	float prev_a; //storage values for incident and reflected waves
	float fs; // sample rate needed for adaptation
	float alpha; //alpha parameter needed for discretization
	float b_coef;
	float a1_coef;
	float a_coef;
};

class IdealVoltageSource : public wdfNode
{
public:
	IdealVoltageSource() : wdfNode("IdealVoltageSource")
	{
		calcImpedance();
	}

	void calcImpedance() {}

	void setVoltage(float V)
	{
		Vs = V;
	}

	void calcIncidentWave(float downWave)
	{
		a = downWave;
	}

	float calcReflectedWave()
	{
		b = 2.0 * Vs - a;
		return b;
	}

private:
	float Vs;
};

class ResistiveVoltageSource : public wdfNode
{
public:
	ResistiveVoltageSource(float R) : wdfNode("ResistiveVoltageSource"), Vs(0.0), R_ser(R)
	{
		calcImpedance();
	}

	void calcImpedance()
	{
		Rp = R_ser;
	}

	void calcIncidentWave(float downSamp)
	{
		a = downSamp;
	}

	float calcReflectedWave()
	{
		b = Vs;
		return b;
	}

	void setVoltage(float V)
	{
		Vs = V;
	}

	void setResistance(float R)
	{
		if (R == R_ser) return;

		R_ser = R;
		propagateImpedance();
	}

private:
	float Vs;
	float R_ser;
};

class IdealCurrentSource : public wdfNode
{
public:
	IdealCurrentSource() : wdfNode("IdealCurrentSource"), Is(0.0)
	{
		calcImpedance();
	}

	~IdealCurrentSource() {}

	void calcImpedance() {}

	void setCurrent(float I)
	{
		Is = I;
	}

	void calcIncidentWave(float downWave)
	{
		a = downWave;
	}

	float calcReflectedWave()
	{
		b = 2.0 * connectedNode->Rp * Is + a;
		return b;
	}

private:
	float Is;
};


class ResistiveCurrentSource : public wdfNode
{
public:
	ResistiveCurrentSource(float R) : wdfNode("ResistiveCurrentSource"), Is(0.0), R_par(R)
	{
		calcImpedance();
	}

	~ResistiveCurrentSource() {}

	void calcImpedance()
	{
		Rp = R_par;
	}

	void setResistance(float R)
	{
		if (R == R_par) return;

		R_par = R;
		propagateImpedance();
	}

	void setCurrent(float I)
	{
		Is = I;
	}

	void calcIncidentWave(float downWave)
	{
		a = downWave;
	}

	float calcReflectedWave()
	{
		b = 2 * Rp * Is;
		return b;
	}

private:
	float Is;
	float R_par;
};

class Switch : public wdfNode
{
public:
	Switch() : wdfNode("Switch"), switchClosed(true) {}
	~Switch() {}

	void calcImpedance() {} //don't calculate impedance, it must be root

	void calcIncidentWave(float downSamp)
	{
		a = downSamp;
	}

	float calcReflectedWave()
	{
		b = switchClosed ? -a : a;
		return b;
	}

	void closeSwitch(bool pos)
	{
		switchClosed = pos;
	}

private:
	bool switchClosed;
};

class TwoPortAdaptor : public wdfNode
{
public:
	TwoPortAdaptor(wdfNode* childPort, std::string type) :
		wdfNode(type), childPort(childPort)
	{
		childPort->connectToNode(this);
	}

	virtual ~TwoPortAdaptor() {}

protected:
	wdfNode* childPort;
};


class Inverter : public TwoPortAdaptor
{
public:
	Inverter(wdfNode* childPort) : TwoPortAdaptor(childPort,"Inverter")
	{
		calcImpedance();
	}

	~Inverter() {}

	void calcImpedance()
	{
		Rp = childPort->Rp;
	}

	void calcIncidentWave(float downWave)
	{
		a = downWave;
		childPort->calcIncidentWave(-a);
	}

	float calcReflectedWave()
	{
		b = -childPort->calcReflectedWave();
		return b;
	}
};

class TwoPortParallel : public TwoPortAdaptor
{
public:
	TwoPortParallel(wdfNode* childPort) : TwoPortAdaptor(childPort, "TwoPortParallel")
	{
		calcImpedance();
	}

	~TwoPortParallel() {}

	void calcImpedance()
	{
		Rp = childPort->Rp;
	}

	void calcIncidentWave(float downWave)
	{
		a = downWave;
		childPort->calcIncidentWave(a);
	}

	float calcReflectedWave()
	{
		b = childPort->calcReflectedWave();
		return b;
	}
};

class IdealTransformer : public TwoPortAdaptor
{
public:
	IdealTransformer(wdfNode* childPort, float turnRatio) : TwoPortAdaptor(childPort,"IdealTransformer"), turnRatio(turnRatio)
	{
		calcImpedance();
	}

	void setTurnRatio(float n)
	{
		turnRatio = n;
		calcImpedance();
	}

	void calcImpedance()
	{
		Rp = turnRatio * turnRatio * childPort->Rp;
	}

	void calcIncidentWave(float downWave)
	{
		a = downWave;
		childPort->calcIncidentWave(a / turnRatio);
	}

	float calcReflectedWave()
	{
		b = childPort->calcReflectedWave() * turnRatio;
		return b;
	}

private:
	float turnRatio;

};

class ActiveTransformer : public TwoPortAdaptor
{
public:
	ActiveTransformer(wdfNode* childPort, float turn1, float turn2) : TwoPortAdaptor(childPort, "ActiveTransformer"), turn1(turn1), turn2(turn2)
	{
		calcImpedance();
	}

	void setTurnRatios(float n1, float n2)
	{
		turn1 = n1;
		turn2 = n2;
		calcImpedance();
	}

	void calcImpedance()
	{
		Rp = turn1 * turn2 * childPort->Rp;
	}

	void calcIncidentWave(float downWave)
	{
		a = downWave;
		childPort->calcIncidentWave(a / turn1);
	}

	float calcReflectedWave()
	{
		b = childPort->calcReflectedWave() * turn1;
		return b;
	}

private:
	float turn1, turn2;

};



class ThreePortAdaptor : public wdfNode
{
public:
	ThreePortAdaptor(wdfNode* leftPort, wdfNode* rightPort, std::string type) :
		wdfNode(type), leftPort(leftPort), rightPort(rightPort)
	{
		leftPort->connectToNode(this);
		rightPort->connectToNode(this);
	}

	virtual ~ThreePortAdaptor() {}

protected:
	wdfNode* leftPort;
	wdfNode* rightPort;
};

class Parallel : public ThreePortAdaptor
{
public:
	Parallel(wdfNode* leftPort, wdfNode* rightPort) : ThreePortAdaptor(leftPort, rightPort, "Parallel")
	{
		calcImpedance();
	}

	~Parallel() {}

	void calcImpedance()
	{
		Rp = (leftPort->Rp * rightPort->Rp) / (leftPort->Rp + rightPort->Rp);
		gammaLeft = Rp / leftPort->Rp;
		gammaRight = Rp / rightPort->Rp;
	}

	void calcIncidentWave(float downWave)
	{
		leftPort->calcIncidentWave(downWave + (rightPort->b - leftPort->b) * gammaRight);
		rightPort->calcIncidentWave(downWave + (rightPort->b - leftPort->b) * -gammaLeft);

		a = downWave;
	}

	float calcReflectedWave()
	{
		b = gammaLeft * leftPort->calcReflectedWave() + gammaRight * rightPort->calcReflectedWave();
		return b;
	}

private:
	float gammaLeft, gammaRight;
};

class Series : public ThreePortAdaptor
{
public:
	Series(wdfNode* leftPort, wdfNode* rightPort) : ThreePortAdaptor(leftPort, rightPort, "Series")
	{
		calcImpedance();
	}

	~Series() {}

	void calcImpedance()
	{
		Rp = leftPort->Rp + rightPort->Rp;
		gammaLeft = leftPort->Rp / Rp;
		gammaRight = rightPort->Rp / Rp;
	}

	void calcIncidentWave(float downWave)
	{
		//port->incident = port->b
		leftPort->calcIncidentWave(leftPort->b - gammaLeft * (downWave + leftPort->b + rightPort->b));
		rightPort->calcIncidentWave(rightPort->b - gammaRight * (downWave + leftPort->b + rightPort->b));

		a = downWave;
	}

	float calcReflectedWave()
	{
		//port->calcreflectedwave = port->a
		b = -1.0 * (leftPort->calcReflectedWave() + rightPort->calcReflectedWave());
		return b;
	}

private:
	float gammaLeft, gammaRight;
};

class RtypeAdaptor : public wdfNode
{
public:
	RtypeAdaptor(int numPorts, wdfNode** downPorts) : wdfNode("R-type Adaptor"), numPorts(numPorts), downPorts(downPorts)
	{
		Rp_down = new float[numPorts+1];
		a_ = new float[numPorts+1];
		b_ = new float[numPorts+1];
		S_matrix = new float* [numPorts+1];
		//double duty for loop:
		for (int i = 0; i < numPorts+1; i++)
		{
			S_matrix[i] = new float[numPorts+1]; //fill out S_matrix initialization
			b_[i] = 0.0f;
			a_[i] = 0.0f;
		}

		for (int i = 0; i < numPorts; i++)
			downPorts[i]->connectToNode(this); //connect downports to R_adaptor
	}

	~RtypeAdaptor()
	{
		delete Rp_down;
		delete a_;
		delete b_;
		for (int i = 0; i < numPorts+1; i++)
			delete[] S_matrix[i];
		delete[] S_matrix;
	}

	void calcImpedance() {}

	void setSMatrixData(float** S_)
	{
		S_matrix = S_;
	}

	float** getSMatrixData()
	{
		return S_matrix;
	}

	void calcIncidentWave(float downWave)
	{
		a_[0] = downWave;
		b_[0] = 0.0f;
		for (int i = 0; i < numPorts; i++)
		{
			a_[i+1] = downPorts[i]->b;
			b_[i+1] = 0.0f;
		}

		RtypeScatter(numPorts+1, S_matrix, a_, b_);

		for (int i = 0; i < numPorts; i++)
		{
			downPorts[i]->calcIncidentWave(b_[i+1]);
		}
		a = a_[0];
	}

	float calcReflectedWave()
	{
		float* S_0 = S_matrix[0];

		for (int i = 0; i < numPorts; i++)
		{
			b += S_0[i+1] * downPorts[i]->calcReflectedWave();
		}

		return b;
	}


protected:
	int numPorts; //number of ports connected to RtypeAdaptor
	float* Rp_down; // array of port resistances
	float** S_matrix; //square matrix representing S
	wdfNode** downPorts; //array of ports connected to RtypeAdaptor
	float* a_; //temp matrix of inputs to Rport
	float* b_; //temp matrix of outputs from Rport
};

class RootRtypeAdaptor : public wdfNode
{
public:
	RootRtypeAdaptor(int numPorts, wdfNode** dP) : wdfNode("Root R-type Adaptor"), numPorts(numPorts), downPorts(dP)
	{
		a_ = new float[numPorts];
		b_ = new float[numPorts];
		S_matrix = new float* [numPorts];
		//double duty for loop:
		for (int i = 0; i < numPorts; i++)
		{
			S_matrix[i] = new float[numPorts]; //fill out S_matrix initialization
			downPorts[i]->connectToNode(this); //connect downports to R_adaptor
			b_[i] = 0.0f;
			a_[i] = 0.0f;
		}
		calcImpedance();
	}

	~RootRtypeAdaptor()
	{
		delete a_;
		delete b_;
		for (int i = 0; i < numPorts; i++)
			delete[] S_matrix[i];
		delete[] S_matrix;
	}

	void calcImpedance() {}

	void setSMatrixData(float** S_)
	{
		S_matrix = S_;
	}

	float** getSMatrixData()
	{
		return S_matrix;
	}

	void calcIncidentWave(float downWave)
	{
		for (int i = 0; i < numPorts; i++)
		{
			a_[i] = downPorts[i]->b;
			b_[i] = 0.0f;
		}

		RtypeScatter(numPorts, S_matrix, a_, b_);

		for (int i = 0; i < numPorts; i++)
		{
			downPorts[i]->calcIncidentWave(b_[i]);
		}
	}

	float calcReflectedWave()
	{
		for (int i = 0; i < numPorts; i++)
		{
			downPorts[i]->calcReflectedWave();
		}

		return 0;
	}
protected:
	int numPorts; //number of ports connected to RtypeAdaptor
	float** S_matrix; //square matrix representing S
	wdfNode** downPorts; //array of ports connected to RtypeAdaptor
	float* a_; //temp matrix of inputs to Rport
	float* b_; //temp matrix of outputs from Rport
};

class Diode : public wdfNode
{
public:
	//Is: reverse saturation current
	//Vt: thermal voltage
	Diode(DiodeModel d) : wdfNode("Diode"), d(d) {}
	~Diode() {}

	void calcImpedance() {}

	void calcIncidentWave(float downWave)
	{
		a = downWave;
	}

	float calcReflectedWave()
	{
		b = a + 2 * connectedNode->Rp * d.Is - 2 * d.Vt * d.nD * omega4(logf((connectedNode->Rp * d.Is) / (d.nD * d.Vt)) + (a + connectedNode->Rp * d.Is) / (d.nD * d.Vt));
		return b;
	}

private:
	const DiodeModel d; // get diode model from wdf_utils
};

class DiodePair : public wdfNode
{
public:
	DiodePair(DiodeModel d) : wdfNode("DiodePair"), d(d) {}
	~DiodePair() {}

	void calcImpedance() {}

	void calcIncidentWave(float downWave)
	{
		a = downWave;
	}

	float calcReflectedWave()
	{
		float lambda = static_cast<float>(signum(a));
		b = a + 2 * lambda * (connectedNode->Rp * d.Is - d.Vt * d.nD * omega4(logf(connectedNode->Rp * d.Is / (d.nD * d.Vt)) + (lambda * a + connectedNode->Rp * d.Is) /(d.nD * d.Vt)));
		return b;
	}
private:
	const DiodeModel d;
};


#endif /* WDF_H_ */
