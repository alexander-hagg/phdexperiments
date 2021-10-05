#include <boost/shared_ptr.hpp>
#include <iostream>
#include <time.h>

#include "fastsim.hpp"
#include "mlp.hpp"

int main(int argc, char* argv[])
{
    using namespace fastsim;

    if (argc < 7) {
        std::cerr << "This program cannot run: you need to provide a XML file (e.g.," << argv[0] << " worlds/example.xml), the number of time steps, a task ID, boolean flag whether to use goal sensors, a number of hidden neurons, debug flag, and a number of weights" << std::endl;
        exit(1);    
    }   

    int fcnParms = 8;
    fastsim::Settings settings(argv[1]);
    boost::shared_ptr<Map> map = settings.map();
    boost::shared_ptr<Robot> robot = settings.robot();
    std::string tmpdir = boost::lexical_cast<std::string>(argv[2]);
    int timesteps = boost::lexical_cast<int>(argv[3]);
    std::string coreID = boost::lexical_cast<std::string>(argv[4]);
    bool useGoalFinders = boost::lexical_cast<bool>(argv[5]);    
    int numHiddenNeurons = boost::lexical_cast<int>(argv[6]);
    int debug = boost::lexical_cast<bool>(argv[7]);
    int numInputs = 3;
    if (useGoalFinders) { 
        numInputs = 7;
    }
    int numLayers = 2;

    using namespace nn;
    Mlp<Neuron<PfWSum<>, AfTanhNoBias<> >, Connection<> > nn(numInputs, numHiddenNeurons, 2);
    std::vector<float> in(numInputs);

    // weights
    std::vector<float> w(nn.get_nb_connections());
    if (debug) std::cout<<"argc - fcnParms :" << argc - fcnParms << " --- nb connections:" << nn.get_nb_connections() << std::endl;
    
    assert(argc - fcnParms == nn.get_nb_connections());    
    for (size_t i = fcnParms; i < argc; ++i)
      w[i - fcnParms] = boost::lexical_cast<float>(argv[i]);
    nn.set_all_weights(w);
    nn.init();
    
    boost::shared_ptr<Display> d;
    if (settings.display())
        d = boost::shared_ptr<Display>(new Display(map, *robot));

    // Trajectory output is saved in filename dependent on the coreID to allow multiple runs
    std::string trajFilename = tmpdir + "/traj.dat" + coreID;
    std::ofstream traj(trajFilename.c_str());

    for (int k = 0; k < timesteps; ++k) {
        if (debug) std::cout << "Time step: " << k << "/" << timesteps << std::endl;
        // lasers
        if (debug) std::cout << "Number of lasers: " << robot->get_lasers().size() << std::endl;
        for (size_t j = 0; j < robot->get_lasers().size(); ++j) {
            float d = robot->get_lasers()[j].get_dist();
            float range = robot->get_lasers()[j].get_range();
            in[j] = (d == -1 ? -1 : 1 - 2 * d / range);
        }
        if (debug) std::cout << "Laser values read" << std::endl;
        if (debug) std::cout << useGoalFinders << std::endl;
        
         // radar    
        if (useGoalFinders) { 
            int numSlices = robot->get_radars()[0].get_nb_slices();
            if (debug) std::cout << "Radar vector size: " << numSlices << std::endl;
            for (size_t j = 0; j < numSlices; ++j) {
                in[robot->get_lasers().size() + j] = -1.0f;
            }
            if (debug) std::cout << "Reading radar values" << std::endl;
            int s = robot->get_radars()[0].get_activated_slice();
            if (debug) std::cout << "Radar values read, activated slice: " << s << std::endl;
            in[robot->get_lasers().size() + s] = 1.0f;
            if (debug) std::cout << "Radar values assigned" << std::endl;
        } else {
            if (debug) std::cout << "Radar values selected to be ignored" << std::endl;
        }
        
        if (debug) std::cout << "Stepping network" << std::endl;
        if (debug) std::cout << "Input: " << in << std::endl;
        for (int j = 0; j < numLayers; ++j){
            nn.step(in);
            if (debug) std::cout << nn.get_outf(0) << std::endl;
        }

        if (debug) std::cout << "Stepping network done" << std::endl;
        
        // move the robot
        float out1 = nn.get_outf(0);
        float out2 = nn.get_outf(1);
        robot->move(out1 * 2, out2 * 2, map);
        if (debug) std::cout << "Moved robot" << std::endl;

        traj << robot->get_pos().x() << " " 
             << robot->get_pos().y() << " " 
             << robot->get_pos().theta() << std::endl;
        // update the display
        if (settings.display())
            d->update();

    }
    return 0;
}
