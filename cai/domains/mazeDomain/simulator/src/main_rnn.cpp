#include <boost/shared_ptr.hpp>
#include <iostream>
#include <time.h>

#include "fastsim.hpp"
#include "mlp.hpp"

int main(int argc, char* argv[])
{
    using namespace fastsim;
    // if (argc != 35 + 2) {
    //     std::cerr << "This program cannot run: you need to provide a XML file (e.g.," << argv[0] << " worlds/example.xml) and 35 weights" << std::endl;
    //    exit(1);
    //  }
    fastsim::Settings settings(argv[1]);
    boost::shared_ptr<Map> map = settings.map();
    boost::shared_ptr<Robot> robot = settings.robot();
    int timesteps = boost::lexical_cast<int>(argv[2]);
    std::string coreID = boost::lexical_cast<std::string>(argv[3]);

    using namespace nn;
    Mlp<Neuron<PfWSum<>, AfTanhNoBias<> >, Connection<> > nn(4 + 12, 3, 2);
    std::vector<float> in(4 + 12);
    // Mlp<Neuron<PfWSum<>, AfTanhNoBias<> >, Connection<> > nn(4, 3, 2);
    // std::vector<float> in(4);

    // weights
    std::vector<float> w(nn.get_nb_connections());
    // std::cout<<"argc -2 :" << argc - 2 << " --- nb connections:" << nn.get_nb_connections() << std::endl;
    
    assert(argc - 4 == nn.get_nb_connections());    
    //srand(time(0));
    for (size_t i = 4; i < argc; ++i)
      w[i - 4] = boost::lexical_cast<float>(argv[i]);
    //rand()/(RAND_MAX+1.0)*4-2.0;
    nn.set_all_weights(w);
    nn.init();
    //std::ofstream nn_ofs("nn.dot");
    //nn.write(nn_ofs);
    
    boost::shared_ptr<Display> d;
    if (settings.display())
        d = boost::shared_ptr<Display>(new Display(map, *robot));

    std::string trajFilename = "traj.dat" + coreID;
    std::ofstream traj(trajFilename.c_str());
    for (int k = 0; k < timesteps; ++k) {
        // lasers
        for (size_t j = 0; j < robot->get_lasers().size(); ++j) {
            float d = robot->get_lasers()[j].get_dist();
            float range = robot->get_lasers()[j].get_range();
            in[j] = (d == -1 ? -1 : 1 - 2 * d / range);
        }
        // radar
        for (size_t j = 0; j < 4; ++j)
            in[robot->get_lasers().size() + j + 1] = -1.0f;
        int s = robot->get_radars()[0].get_activated_slice();
        
        in[robot->get_lasers().size() + 1] = 1.0f;

        for (int j = 0; j < 3; ++j)
            nn.step(in);

        // move the robot
        float out1 = nn.get_outf(0);
        float out2 = nn.get_outf(1);
        robot->move(out1 * 2, out2 * 2, map);
        //std::cout<<"out:"<<out1<<" "<<out2<<std::endl;

        traj << robot->get_pos().x() << " " 
             << robot->get_pos().y() << " " 
             << robot->get_pos().theta() << std::endl;
        // update the display
        if (settings.display())
            d->update();

    }
    return 0;
}
