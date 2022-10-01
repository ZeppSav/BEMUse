#ifndef SOLVER_SETUP_H
#define SOLVER_SETUP_H

#include <QWidget>

//--- BEMUse Objects
#include "src/Solver/Aerodynamic_Solver.h"
#include "src/Solver/Hydrodynamic_Solver.h"

namespace Ui {
class Solver_Setup;
}

class Solver_Setup : public QWidget
{
    Q_OBJECT

    //--- Solver object
    BEMUse::Solver *Solver;

public:
    explicit Solver_Setup(QWidget *parent = nullptr);
    ~Solver_Setup();

signals:
    void Pass_Solver(BEMUse::Solver *S, StdVector Freqs, StdVector Betas);

private slots:

    void on_Conf_Solver_clicked();
//    void on_IRFToggle_toggled(bool checked)     {IRF = checked;}

private:
    Ui::Solver_Setup *ui;
};

#endif // SOLVER_SETUP_H
