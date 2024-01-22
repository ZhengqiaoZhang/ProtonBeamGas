#include "HepMC3/GenEvent.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/ReaderAsciiHepMC2.h"


int transFormat()
{

HepMC3::ReaderAsciiHepMC2 inputA("output.hepmc"); 

if(inputA.failed()) return 1;

HepMC3::WriterAscii outputA("sample1.Asciiv3.hepmc"); 

if(outputA.failed()) return 2;

while( !inputA.failed() )
{
HepMC3::GenEvent evt(HepMC3::Units::GEV,HepMC3::Units::MM); 
inputA.read_event(evt);

if( inputA.failed() ) {printf("End of file reached. Exit.\n"); break;} 

outputA.write_event(evt);
evt.clear();
          }
inputA.close();
outputA.close(); 
return 0;

}
