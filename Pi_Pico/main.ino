extern "C" 
{ 
#include "Subsystem.h" 
}

int i=0;

void setup()
{
    Serial.begin(9600);
    Subsystem_initialize();
}

void loop()
{
    if (i > 1.06/0.001)
    {
        Subsystem_U.G_fuel = 0.0023;
    }
    if (i > 5/0.001)
    {
        Subsystem_U.G_fuel = 0.0045;
    }
    if (i > 8/0.001)
    {
        Subsystem_U.G_fuel = 0.0085;
    }
    Subsystem_step();
    Serial.print(Subsystem_Y.P);
    Serial.print('\t');
    Serial.print(i*0.001);
    Serial.println();
    i++;
}
