#include "config.hpp"
#include "elastic_wave3D.hpp"
#include "mfem.hpp"
#include "parameters.hpp"
#include "utilities.hpp"

using namespace std;
using namespace mfem;



int main(int argc, char *argv[])
{
  if (argc == 1) // no arguments
  {
    cout << "\nGet help with:\n" << argv[0] << " -h\n" << endl;
    return 0;
  }

#if defined(DEBUG_WAVE)
  cout << "****************************\n";
  cout << "*     DEBUG VERSION        *\n";
  cout << "****************************\n";
#endif

  try
  {
    StopWatch chrono;
    chrono.Start();

    Parameters param;
    param.init(argc, argv);

//    show_SRM_damp_weights(param);

    ElasticWave2D elwave(param);
    elwave.run();

    cout << "\nTOTAL TIME " << chrono.RealTime() << " sec\n" << endl;
  }
  catch (int ierr)
  {
    return ierr;
  }
  catch (...)
  {
    cerr << "\nEXCEPTION\n";
  }

  return 0;
}
