#ifndef GRID_OPTIONS_H
#define GRID_OPTIONS_H

#include <QWidget>
#include <QList>
#include <QLabel>
#include <QDoubleSpinBox>
#include <QSpinBox>
#include <QKeyEvent>
#include <QCheckBox>
#include <QPushButton>

// Include geometry options
#include "src/Boundary/Ellipsoid.h"
#include "src/Boundary/Volume_of_Revolution.h"
#include "src/Boundary/Horizontal_Volume_of_Revolution.h"
#include "src/Boundary/Triple_Spar.h"
#include "src/Boundary/Ship_Hulls.h"
#include "src/Boundary/STL_Geo.h"
#include "src/Boundary/GDF_Geo.h"
#include "src/Boundary/MAR_Geo.h"
#include "src/Boundary/PNL_Geo.h"
#include "src/Boundary/FOWT_Platforms.h"
#include "src/Boundary/Barge.h"
#include "src/Boundary/Thin_Disc.h"

#include "src/Boundary/Wing.h"

enum Geo    {   No_Geo,
                Semi_Ellipsoid,
                Ellipsoid,
                Half_Cylinder,
                Barge,
                Tapered_SparBuoy,
                OC3_SparBuoy,
                Triple_Spar,
                SemiSub_OC4,
                Wigley_Hull};

namespace Ui {
class Grid_Options;
}

static int NDimLabels = 10;
static int NDiscLabels = 10;

class Grid_Options : public QWidget
{
    Q_OBJECT

    // Boundary options
    Geo GeoType = No_Geo;

    // Boundary object
    BEMUse::Boundary *Boundary = NULL;

    // Visualisation lists.
    std::vector<QString>  DimLabel_Text, DiscLabel_Text;
//    std::vector<bool>     Toggle_Dim, Toggle_Disc;
    std::vector<bool>     Toggle_Dim, Toggle_DimAux, Toggle_DimExt;
    std::vector<bool>     Toggle_Disc, Toggle_DiscAux, Toggle_DiscExt;

//    QList<QLabel*>          DimLabList, DiscLabList;
//    QList<QDoubleSpinBox*>  DimSpinList;
    QList<QSpinBox*>        DiscSpinList;

    // New formulation

    int NDims = 0, NDimsAux=0, NDimsExt = 0;
    int NDisc = 0, NDiscAux=0, NDiscExt = 0;

//    std::vector<Real>       Dims, DimsAux, DimsExt;
//    std::vector<int>        Disc, DiscAux, DiscExt;

    QList<QLabel*>          DimLabs, DimLabsAux, DimLabsExt;
    QList<QLabel*>          DiscLabs, DiscLabsAux, DiscLabsExt;

    QList<QDoubleSpinBox*>  DimSpins, AuxDimSpins, ExtDimSpins;
    QList<QSpinBox*>        DiscSpins, AuxDiscSpins, ExtDiscSpins;
    QList<QPushButton*>     Buttons;

    QList<QCheckBox*>       CheckBoxes;

    void Hide_Labels();
    void Hide_Buttons();

    std::vector<Real> DimInitVals;
    std::vector<int> DiscInitVals;

    //--- Visualisations
//    void Lighten_Button(QPushButton *Button);
    void Darken_Button(QPushButton *Button);

public:
    explicit Grid_Options(QWidget *parent = nullptr);
    ~Grid_Options()     {delete ui;}

signals:
    void Pass_Geometry(BEMUse::Boundary *B);

protected:

    //--- Virtual QOpenGLWidget functions
//    void Reset_Labels();
//    void Toggle_Dimension_Inputs();
//    void Toggle_Discretisation_Inputs();

private slots:
    void on_ButtonSemiEllipsoid_clicked();
    void on_ButtonEllipsoid_clicked();
    void on_ButtonCylinder_clicked();
    void on_ButtonBarge_clicked();
    void on_ButtonSparBuoy_clicked();
    void on_ButtonOC3_clicked();
    void on_ButtonTripleSpar_clicked();
    void on_ButtonOC4_clicked();
    void on_ButtonWigley_clicked();

    void on_GenerateButton_clicked();

private:
    Ui::Grid_Options *ui;
};

#endif // GRID_OPTIONS_H
