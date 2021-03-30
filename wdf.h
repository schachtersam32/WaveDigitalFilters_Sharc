#include <matrix.h>
#include <vector>
#include <stddef.h>
#include <stdlib.h>

class wdfTree //a special instance of this class is created for each circuit to be run
{
public:
	wdfTree();
	virtual ~wdfTree() {}

	virtual void setInputVoltage(float inSamp) = 0; //user defines where input voltage is sent
	virtual float getOutputVoltage() = 0; //user defines where to take voltage from
	
	float processSample() //processes tree sample by sample
	{
		root->calcIncidentWave(bottomAdaptor->calcReflectedWave()); //calculate waves incident to root by reflecting from leaves
		bottomAdaptor->calcIncidentWave(root->calcReflectedWave()); //calculate waves reflected from root incident to leaves
	}

	float setSampleRate(float SR) { sr = SR; } //set system samplerate
	float setAlpha(float a) { alpha = a; } //set circuit alpha parameter
protected:
	wdfNode* root;
	wdfNode* bottomAdaptor;
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
	wdfNode(std::string type) : wdfPort(type) {} //initialize a node with a port
	virtual ~wdfNode() {}

	void connectToNode(wdfPort* node) { connectedNode = node; } //connect ports together

	void propagateImpedance() override
	{
		calcImpedance();

		if (connectedNode != nullptr)
			connectedNode->propagateImpedance();
	}

protected:
	wdfPort* connectedNode = nullptr;

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
	R_val = 0.1;
};

class Capacitor : public wdfNode
{
public:
	Capacitor(float C, float fs, float alpha) : wdfNode("Capacitor"),
		C_val(C), prev_a(0.0), prev_b(0.0), alpha(alpha), fs(fs)
	{
		calcImpedance();
	}

	~Capacitor() {}

	void setCapacitance(float C)
	{
		if (C = C_val)
			return;

		C_val = C;
		propagateImpedance();
	}

	void calcImpedance() override
	{
		Rp = 1 / (C_val * fs * (1 + alpha));
	}

	void calcIncidentWave(float downSamp)
	{
		a = downWave;
		prev_a = a;
	}

	float calcReflectedWave()
	{
		b = b_coef * prev_b + a_coef * prev_a;
		prev_b = b;
		return b;
	}

	float b_coef = (1 - alpha) / 2.0;
	float a_coef = (1 + alpha) / 2.0;
	
private:
	float C_val = 0.0;
	float prev_a, prev_b; //storage values for incident and reflected waves
	float fs; // sample rate needed for adaptation
	float alpha; //alpha parameter needed for discretization
};

class Inductor : public wdfNode
{
public:
	Inductor(float L, float fs, float alpha) : wdfNode("Inductor"),
		L_val(L), prev_a(0.0), prev_b(0.0), alpha(alpha), fs(fs)
	{
		calcImpedance();
	}

	~Inductor() {}

	void setInductance(float L)
	{
		if (L = L_val)
			return;

		L_val = L;
		propagateImpedance();
	}

	void calcImpedance() override
	{
		Rp = (L_val * fs * (1 + alpha));
	}

	void calcIncidentWave(float downSamp)
	{
		a = downWave;
		prev_a = a;
	}

	float calcReflectedWave()
	{
		b = b_coef * prev_b - a_coef * prev_a;
		prev_b = b;
		return b;
	}

	float b_coef = (1 - alpha) / 2.0;
	float a_coef = (1 + alpha) / 2.0;

private:
	float L_val = 0.0;
	float prev_a, prev_b; //storage values for incident and reflected waves
	float fs; // sample rate needed for adaptation
	float alpha; //alpha parameter needed for discretization
};

class IdealVoltageSource : public wdfNode
{
public:
	IdealVoltageSource() : wdfNode("IdealVoltageSource")
	{
		calcImpedance();
	}

	void calcImpedance() {}

	void setVoltage (float V)
	{
		Vs = V;
	}

	void calcIncidentWave(float downWave) override
	{
		a = downWave;
	}

	void calcReflectedWave() override
	{
		b = 2 * Vs - a;
		return b;
	}

private:
	float Vs;
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
	}

	float calcReflectedWave()
	{
		b_T = -1.0 * gamma * (leftPort->a - a);
		b = b_t + leftPort->reflected();
		return b;
	}

	void calcIncidentWave(float downWave)
	{
		b_T = -1.0 * gamma * (leftPort->a - a);
		leftport->incident(b_T + downWave);
		rightPort->incident(leftPort->b + rightPort->a - rightPort->a);

		a = downWave;
	}

private:
	float b_T = 0.0; //scattering parameter based on Fettweis
	float gamma = 0.0;
};

class Series : wdfThreePortAdaptor
{
public:
	Series(wdfNode* leftPort, wdfNode* rightPort) : wdfThreePortAdaptor(leftPort, rightPort, "Series")
	{
		calcImpedance();
	}

	~Series() {}

	void calcImpedance()
	{
		R_left = leftPort->Rp;
		R_right = rightPort->Rp;

		Rp = R_left + R_right;
		gamma = R_left / Rp;
	}

	float calcReflectedWave()
	{
		b = -1.0 * (leftPort->reflected() + rightPort->reflected());
		return b;
	}

	void calcIncidentWave(float downWave)
	{
		leftPort->incident(rightPort->a - gamma * (leftPort->a + rightPort->a + downWave));
		rightPort->incident(-1.0 * (rightPort->b + downWave));

		a = downWave;
	}
private:
	float gamma = 0.0;
};