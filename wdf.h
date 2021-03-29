#include <matrix.h>
#include <vector>
#include <stddef.h>
#include <stdlib.h>

class wdfTree //a special instance of this class is created for each circuit to be run
{
public:
	wdfTree() {}
	virtual ~wdfTree() {}

	virtual float processSample(float inSamp) = 0;
	float setSampleRate(float SR) { sr = SR; }
	float setAlpha(float a) { alpha = a; }
private:
	float sr;
	float alpha;
};
class wdfPort
{
protected:
	float a; //incident wave
	float b; //reflected wave
	float Rp; //port impedance
public:
	wdfPort(std::string portType) : portType(portType) {}
	virtual ~wdfPort() {}

	virtual void calculateImpedance() = 0; //calculates the port impedance
	virtual void propagateImpedance() = 0; //propagate port impedance up the tree

	virtual void calculateIncidentWave(float x) = 0; //find the wave leaving the port
	virtual void calculateReflectedWave() = 0; //find the wave entering the port

	float getPortVoltage()
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
	wdfNode(std::string type) : wdfPort(type) {}
	virtual ~wdfNode() {}

	void connectToNode(wdfPort* node) { connectedNode = node; } //connect ports together

	void propagateImpedance() override
	{
		calcImpedance();

		if (connectedNode != nullptr)
			connectedNode->propagateImpedance();
	}

protected:
	wdfPort* connectedNode;

};

class Resistor : public wdfNode
{
public:
	Resistor(float R) : wdfNode("Resistor"), R_val(R)
	{
		calcImpedance();
	}

	~Resistor() {}

private:
	R_val = 0.1;
};