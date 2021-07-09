// cmptl
// Created by parnet on 08.06.19.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_BRAIDUSAGETIMER_H
#define UG_PLUGIN_XBRAIDFORUG4_BRAIDUSAGETIMER_H

#include <chrono>
#include <ctime>
#include <ratio>
#include <vector>


class BraidUsageTimer {
public: // todo set better modes
    std::chrono::high_resolution_clock::time_point t0;
    std::chrono::high_resolution_clock::time_point t1;

    double time = 0;
    int usage = 0;

    BraidUsageTimer() = default;

    void start() {
        t0 = std::chrono::high_resolution_clock::now();
    }

    void stop() {
        t1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> difference =
                std::chrono::duration_cast < std::chrono::duration < double >> (t1 - t0);
        time += difference.count();
        usage++;
    }

    double getTime() const {
        return time;
    }

    double getUsage() const {
        return usage;
    }

    double getAverageTime() const{
        return time / usage;
    }
};

enum Observer {
    T_INIT = 0, T_CLONE, T_FREE, T_ACCESS, T_SUM, T_SEND, T_RECV, T_STEP,
    T_RESIDUAL, T_NORM, T_COARSEN, T_REFINE, T_SYNC, N_OBSERVER
};

std::string ObserverNames[] = {"init", "clone", "free", "access", "sum", "send", "recv", "step",
                               "residual", "norm", "coarsen", "refine", "sync", "counter"};

enum LevelObserver {
    TL_STEP = 0,
    TL_RESIDUAL,
    TL_ASSEMBLE_OP,
    TL_ASSEMBLE_RHS,
    TL_ADAPTIVE_TOL,
    TL_SOLVE,
    NL_LEVELOBSERVER
};
int cLevelObserver = 6;
std::string LevelObserverNames[] = {"step", "residual", "assemble_op", "assemble_rhs", "adaptive_tol", "solve"};

class BraidTimeLogManager {
    typedef  std::vector<BraidUsageTimer> TimerList;
    typedef std::vector<std::vector<BraidUsageTimer>> TimerMatrix;
public:

    TimerList timer = TimerList();
    TimerMatrix leveltimer = TimerMatrix();

    BraidTimeLogManager() = default;

    explicit BraidTimeLogManager(int maxlevel) {
        for (int i = 0; i < Observer::N_OBSERVER; i++) {
            timer.emplace_back(BraidUsageTimer());
        }


        leveltimer = TimerMatrix(cLevelObserver,TimerList(maxlevel, BraidUsageTimer()));
    }

    BraidUsageTimer &get(Observer o) {
        //std::cout <<"Observer: " <<  o << std::endl;
        return this->timer[o];
    }

    BraidUsageTimer &get(LevelObserver o, int level) {
        //std::cout <<"LObserver: " <<  o << "&" << level << "be"<< std::endl;
        BraidUsageTimer &v = this->leveltimer[o][level];
        //std::cout <<"LObserver: " <<  o << "&" << level << "af"<< std::endl;
        return v;
    };

};

#endif //UG_PLUGIN_XBRAIDFORUG4_BRAIDUSAGETIMER_H
