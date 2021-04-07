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
#include <stddef.h>
#include <stdlib.h>

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
	wdfTree() : alpha(1.0), sr(0.0) {}
	virtual ~wdfTree() {}

	virtual float processSample(float inSamp) = 0; //processes tree sample by sample

	void setSampleRate(float SR) { sr = SR; } //set system samplerate
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

private:
	float R_val;
};

class Capacitor : public wdfNode
{
public:
	Capacitor(float C, float fs, float alpha) : wdfNode("Capacitor"),
		C_val(C), prev_a(0.0), prev_b(0.0), alpha(alpha), fs(fs), b_coef((1 - alpha) / 2.0), a_coef((1 + alpha) / 2.0)
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

private:
	float C_val;
	float prev_a, prev_b; //storage values for incident and reflected waves
	float fs; // sample rate needed for adaptation
	float alpha; //alpha parameter needed for discretization
	const float b_coef;
	const float a_coef;
};

class Inductor : public wdfNode
{
public:
	Inductor(float L, float fs, float alpha) : wdfNode("Inductor"),
		L_val(L), prev_a(0.0), prev_b(0.0), alpha(alpha), fs(fs), b_coef((1 - alpha) / 2.0), a_coef((1 + alpha) / 2.0)
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


private:
	float L_val;
	float prev_a, prev_b; //storage values for incident and reflected waves
	float fs; // sample rate needed for adaptation
	float alpha; //alpha parameter needed for discretization
	const float b_coef;
	const float a_coef;
};

class IdealVoltageSource : public wdfNode
{
public:
	IdealVoltageSource(float Vs) : wdfNode("IdealVoltageSource")
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
		b = 2 * Vs - a;
		return b;
	}

private:
	float Vs;
};

class ResistiveVoltageSource : public wdfNode
{
public:
	ResistiveVoltageSource(float R) : wdfNode("ResistiveVoltageSource"), Vs(0.0), R_ser(1.0) {}

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

class Inverter : public wdfNode
{
public:
	Inverter(wdfNode* port) : wdfNode("Inverter"), port(port)
	{
		port->connectToNode(this);
		calcImpedance();
	}

	~Inverter() {}

	void calcImpedance()
	{
		Rp = port->Rp;
	}

	void calcIncidentWave(float downWave)
	{
		a = downWave;
		port->calcIncidentWave(-a);
	}

	float calcReflectedWave()
	{
		b = -port->calcReflectedWave();
		return b;
	}

private:
	wdfNode* port;
};

class wdfThreePortAdaptor : public wdfNode
{
public:
	wdfThreePortAdaptor(wdfNode* leftPort, wdfNode* rightPort, std::string type) :
		wdfNode(type), leftPort(leftPort), rightPort(rightPort)
	{
		leftPort->connectToNode(this);
		rightPort->connectToNode(this);
	}

	virtual ~wdfThreePortAdaptor() {}

protected:
	wdfNode* leftPort;
	wdfNode* rightPort;
};

class Parallel : public wdfThreePortAdaptor
{
	Parallel(wdfNode* leftPort, wdfNode* rightPort) : wdfThreePortAdaptor(leftPort, rightPort, "Parallel")
	{
		calcImpedance();
	}

	~Parallel() {}

	void calcImpedance()
	{
		const float R_left = leftPort->Rp;
		const float R_right = rightPort->Rp;

		Rp = (R_left * R_right) / (R_left + R_right);
		gamma = Rp / R_right;
		b_T = -1.0 * gamma * (leftPort->a - a);
	}

	float calcReflectedWave()
	{
		b = b_T + leftPort->calcReflectedWave();
		return b;
	}

	void calcIncidentWave(float downWave)
	{
		float left_b = b_T + downWave;
		leftPort->calcIncidentWave(left_b);
		rightPort->calcIncidentWave(left_b + leftPort->a - rightPort->a);

		a = downWave;
	}

private:
	float b_T; //scattering parameter based on Fettweis
	float gamma;
};

class Series : public wdfThreePortAdaptor
{
public:
	Series(wdfNode* leftPort, wdfNode* rightPort) : wdfThreePortAdaptor(leftPort, rightPort, "Series")
	{
		calcImpedance();
	}

	~Series() {}

	void calcImpedance()
	{
		const float R_left = leftPort->Rp;
		const float R_right = rightPort->Rp;

		Rp = R_left + R_right;
		gamma = R_right / Rp;
	}

	float calcReflectedWave()
	{
		b = -1.0 * (leftPort->calcReflectedWave() + rightPort->calcReflectedWave());
		return b;
	}

	void calcIncidentWave(float downWave)
	{
		float right_b = rightPort->a - gamma * (leftPort->a + rightPort->a + downWave);
		rightPort->calcIncidentWave(right_b);
		leftPort->calcIncidentWave(-1.0 * (right_b + downWave));

		a = downWave;
	}
private:
	float gamma;
};


#endif /* WDF_H_ */
