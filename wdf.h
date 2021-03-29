#include <matrix.h>
#include <vector>
#include <stddef.h>
#include <stdlib.h>

class wdfPort
{
protected:
	float a; //incident wave
	float b; //reflected wave
	float Rp; //port impedance
public:
	wdfPort(std::string portType) : portType(portType) {}
	virtual ~wdfPort() {}

	virtual void calculateImpedance() = 0;
	virtual void propagateImpedance() = 0;

	virtual void calculateIncidentWave(float x) = 0;
	virtual void calculateReflectedWave() = 0;

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