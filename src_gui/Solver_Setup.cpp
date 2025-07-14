#include "Solver_Setup.h"
#include "ui_solver_setup.h"
#include <QFileDialog>

Solver_Setup::Solver_Setup(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Solver_Setup)
{
    ui->setupUi(this);

    // Prepare values (cannot currently handle finite depth cases, so deactivate)
    ui->FinDepthButton->hide();
    ui->FinDepthSet->hide();

    // Initial values are set in the GUI for simplicity
}

Solver_Setup::~Solver_Setup()
{
    delete ui;
}

inline Real PiCosFac(int i, int itot)      {return 0.5*(1.0-cos(i*PI/(itot-1)));}
inline Real HalfPiCosFac(int i, int itot)  {return -cos(0.5*PI+0.5*i*PI/(itot-1));}

void Solver_Setup::on_Conf_Solver_clicked()
{
    // The solver parameters have been specified.
    // These are passed to the generated solver object and this is then passed to the main window.
    // The task of the configure configure window is then complete!

    //--- Sim parameters
    std::vector<BEMUse::Parameter> Params;

    //--- Sim flags
    if (ui->IRFToggle->isChecked())     Params.push_back(BEMUse::Parameter("IFR",true));

    //--- Frequency list
    Real F1 = ui->FLowSet->value();
    Real F2 = ui->FHighSet->value();
    int FD = ui->FDiscSet->value();

    StdVector Freqs;
    if (FD==1) Freqs.push_back(F1);
    else if (ui->FDiscLinButton->isChecked())  {for (int i=0; i<FD; i++)    Freqs.push_back(F1 + (F2-F1)*i/(FD-1)); }
    else if (ui->FDiscCos1Button->isChecked()) {for (int i=0; i<FD; i++)    Freqs.push_back(F2 + (F1-F2)*HalfPiCosFac(FD-1-i,FD)); }
    else if (ui->FDiscCos2Button->isChecked()) {for (int i=0; i<FD; i++)    Freqs.push_back(F1 + (F2-F1)*HalfPiCosFac(i,FD)); }

//    for (int i=0; i<Freqs.size(); i++) Freqs[i] *= TwoPIinv;        // Convert

    //--- Wave angle list
    Real B1 = ui->BetaLowSet->value();
    Real B2 = ui->BetaHighSet->value();
    int BD = ui->BetaDiscSet->value();

    StdVector Betas;
    if (BD==1) Betas.push_back(B1);
    else    {
        for (int i=0; i<BD; i++)    Betas.push_back(B1 + (B2-B1)*i/(BD-1));
    }
    for (int i=0; i<Betas.size(); i++) Betas[i] *= D2R;        // Convert
    for (int i=0; i<Betas.size(); i++) {Params.push_back(BEMUse::Parameter("WaveAngle",Real(Betas[i])));}


    //--- Kochin angles
    Params.push_back(BEMUse::Parameter("NKochin",ui->KochinDiscSet->value()));

    //--- Environmental variables
    Params.push_back(BEMUse::Parameter("Density",Real(ui->DensitySpin->value())));
    Params.push_back(BEMUse::Parameter("Gravity",Real(ui->GravitySpin->value())));

    //--- Output names
    std::string OutputName = ui->FileNameEdit->text().toStdString();

    //--- Create solver object
    Solver = new BEMUse::Hydrodynamic_Radiation_Solver();

    //--- Set solver parameters
    Solver->Set_Parameters(Params);
    Solver->Set_OutputFilePath(OutputName);

    //--- Pass this object to the MainWindown
    Pass_Solver(Solver,Freqs,Betas);

    //--- Close window
    this->close();
}
