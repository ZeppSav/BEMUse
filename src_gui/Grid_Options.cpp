#include "Grid_Options.h"
#include "ui_grid_options.h"

Grid_Options::Grid_Options(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Grid_Options)
{
    ui->setupUi(this);

    // Collect all labels, spin boxes and textEdits to allow for simple treatment later
    DimLabs.push_back(ui->DimLabel_1);
    DimLabs.push_back(ui->DimLabel_2);
    DimLabs.push_back(ui->DimLabel_3);
    DimLabs.push_back(ui->DimLabel_4);
    DimLabs.push_back(ui->DimLabel_5);
    DimLabsAux.push_back(ui->DimLabel_6);
    DimLabsAux.push_back(ui->DimLabel_7);
    DimLabsExt.push_back(ui->DimLabel_8);
    DimLabsExt.push_back(ui->DimLabel_9);
    DimLabsExt.push_back(ui->DimLabel_10);

    DimSpins.push_back(ui->DimSpin1);
    DimSpins.push_back(ui->DimSpin2);
    DimSpins.push_back(ui->DimSpin3);
    DimSpins.push_back(ui->DimSpin4);
    DimSpins.push_back(ui->DimSpin5);
    AuxDimSpins.push_back(ui->DimSpin6);
    AuxDimSpins.push_back(ui->DimSpin7);
    ExtDimSpins.push_back(ui->DimSpin8);
    ExtDimSpins.push_back(ui->DimSpin9);
    ExtDimSpins.push_back(ui->DimSpin10);

    DiscLabs.push_back(ui->DiscLabel_1);
    DiscLabs.push_back(ui->DiscLabel_2);
    DiscLabs.push_back(ui->DiscLabel_3);
    DiscLabs.push_back(ui->DiscLabel_4);
    DiscLabs.push_back(ui->DiscLabel_5);
    DiscLabsAux.push_back(ui->DiscLabel_6);
    DiscLabsAux.push_back(ui->DiscLabel_7);
    DiscLabsExt.push_back(ui->DiscLabel_8);
    DiscLabsExt.push_back(ui->DiscLabel_9);
    DiscLabsExt.push_back(ui->DiscLabel_10);

    DiscSpins.push_back(ui->spinBox_1);
    DiscSpins.push_back(ui->spinBox_2);
    DiscSpins.push_back(ui->spinBox_3);
    DiscSpins.push_back(ui->spinBox_4);
    DiscSpins.push_back(ui->spinBox_5);
    AuxDiscSpins.push_back(ui->spinBox_6);
    AuxDiscSpins.push_back(ui->spinBox_7);
    ExtDiscSpins.push_back(ui->spinBox_8);
    ExtDiscSpins.push_back(ui->spinBox_9);
    ExtDiscSpins.push_back(ui->spinBox_10);

    Buttons.push_back(ui->ButtonBarge);
    Buttons.push_back(ui->ButtonCylinder);
    Buttons.push_back(ui->ButtonEllipsoid);
    Buttons.push_back(ui->ButtonSemiEllipsoid);
    Buttons.push_back(ui->ButtonOC3);
    Buttons.push_back(ui->ButtonOC4);
    Buttons.push_back(ui->ButtonSparBuoy);
    Buttons.push_back(ui->ButtonTripleSpar);
    Buttons.push_back(ui->ButtonWigley);

    CheckBoxes.push_back(ui->Flag1);
    CheckBoxes.push_back(ui->Flag2);

    Hide_Labels();
    Hide_Buttons();
}

void Grid_Options::on_GenerateButton_clicked()
{
    // The "generate" button has been clicked.
    // This means we should proceed with the generation and preparation of the geometry
    // Finally the geometry is passed to the BEMUser Interface QObject for the next steps

    switch (GeoType)
    {
    //--- Template geometries
    case No_Geo:            {return;   break;}  // Jump out if no geometry type has been selected.
    case Semi_Ellipsoid:    {Generate_SemiEllipsoid();          break;}
    case Ellipsoid:         {Generate_Ellipsoid();              break;}
    case Half_Cylinder:     {Generate_Cylinder();               break;}
    case Barge:             {Generate_Barge();                  break;}
        // case NWT:               {Generate_Numerical_Wave_Tank();    break;}

        //--- FOWT platforms
    case Tapered_SparBuoy:  {Generate_SparBuoy();               break;}
    case OC3_SparBuoy:      {Generate_OC3SparBuoy();            break;}
    case Triple_Spar:       {Generate_TripleSpar();             break;}
    case SemiSub_OC4:       {Generate_OC4TripleSpar();          break;}

        //--- Ship geometries
    case Wigley_Hull:       {Generate_WigleyHull();             break;}

    default: break;
    }

    // Now setup geometry and close geo dialogue
    Boundary->Set_Parameters(Parameters);
    Boundary->Setup();
    Pass_Geometry(Boundary);
    this->close();
}

void Grid_Options::Hide_Labels()
{
    // Simply hides the labels in case the geometry options has been modified

    for (QLabel* L : DimLabs)               L->hide();
    for (QLabel* L : DimLabsAux)            L->hide();
    for (QLabel* L : DimLabsExt)            L->hide();
    for (QLabel* L : DiscLabs)              L->hide();
    for (QLabel* L : DiscLabsAux)           L->hide();
    for (QLabel* L : DiscLabsExt)           L->hide();

    for (QDoubleSpinBox* L : DimSpins)      L->hide();
    for (QDoubleSpinBox* L : AuxDimSpins)   L->hide();
    for (QDoubleSpinBox* L : ExtDimSpins)   L->hide();
    for (QSpinBox* L : DiscSpins)           L->hide();
    for (QSpinBox* L : AuxDiscSpins)        L->hide();
    for (QSpinBox* L : ExtDiscSpins)        L->hide();

    for (QCheckBox* L : CheckBoxes)         L->hide();

    NDims = 0;
    NDimsAux = 0;
    NDimsExt = 0;

    NDisc = 0;
    NDiscAux = 0;
    NDiscExt = 0;
}

void Grid_Options::Hide_Buttons()
{
    for (QPushButton* B : Buttons)
    {
        QPalette pal = B->palette();
        pal.setColor(QPalette::Button, QColor(Qt::transparent));
        B->setAutoFillBackground(true);
        B->setPalette(pal);
        B->update();
    }
}

void Grid_Options::Darken_Button(QPushButton *Button)
{
    QPalette pal = Button->palette();
    pal.setColor(QPalette::Button, QColor(Qt::gray));
    Button->setAutoFillBackground(true);
    Button->setPalette(pal);
    Button->update();

}

void Grid_Options::on_ButtonSemiEllipsoid_clicked()
{
    // Hemisphere type geometry
    GeoType = Semi_Ellipsoid;
    Hide_Labels();
    Hide_Buttons();
    Darken_Button(ui->ButtonSemiEllipsoid);

    DimLabs[0]->setText("Semiaxis a"); DimLabs[0]->show(); DimSpins[0]->setValue(1.0); DimSpins[0]->show(); NDims++;
    DimLabs[1]->setText("Semiaxis b"); DimLabs[1]->show(); DimSpins[1]->setValue(1.0); DimSpins[1]->show(); NDims++;
    DimLabs[2]->setText("Semiaxis c"); DimLabs[2]->show(); DimSpins[2]->setValue(1.0); DimSpins[2]->show(); NDims++;
    DimLabsExt[0]->setText("Radial factor"); DimLabsExt[0]->show(); ExtDimSpins[0]->setValue(8.0); ExtDimSpins[0]->show(); NDimsExt++;

    DiscLabs[0]->setText("Azimuthal");      DiscLabs[0]->show();    DiscSpins[0]->setValue(32); DiscSpins[0]->show();   NDisc++;
    DiscLabs[1]->setText("Axial");          DiscLabs[1]->show();    DiscSpins[1]->setValue(32); DiscSpins[1]->show();  NDisc++;
    DiscLabsAux[0]->setText("Radial");      DiscLabsAux[0]->show(); AuxDiscSpins[0]->setValue(32); AuxDiscSpins[0]->show();    NDiscAux++;
    DiscLabsExt[0]->setText("Azimuthal");   DiscLabsExt[0]->show(); ExtDiscSpins[0]->setValue(32); ExtDiscSpins[0]->show(); NDiscExt++;
    DiscLabsExt[1]->setText("Axial");       DiscLabsExt[1]->show(); ExtDiscSpins[1]->setValue(32); ExtDiscSpins[1]->show(); NDiscExt++;

    update();
}

void Grid_Options::on_ButtonEllipsoid_clicked()
{
    // Ellipsoid type geometry
    GeoType = Ellipsoid;
    Hide_Labels();
    Hide_Buttons();
    Darken_Button(ui->ButtonEllipsoid);

    DimLabs[0]->setText("Semiaxis a"); DimLabs[0]->show(); DimSpins[0]->setValue(1.0); DimSpins[0]->show(); NDims++;
    DimLabs[1]->setText("Semiaxis b"); DimLabs[1]->show(); DimSpins[1]->setValue(1.0); DimSpins[1]->show(); NDims++;
    DimLabs[2]->setText("Semiaxis c"); DimLabs[2]->show(); DimSpins[2]->setValue(1.0); DimSpins[2]->show(); NDims++;
    DimLabs[3]->setText("Depth");      DimLabs[3]->show(); DimSpins[3]->setValue(0.0); DimSpins[3]->show(); NDims++;
    DimLabsExt[0]->setText("Radius");  DimLabsExt[0]->show(); ExtDimSpins[0]->setValue(15.0); ExtDimSpins[0]->show(); NDimsExt++;

    DiscLabs[0]->setText("Azimuthal");      DiscLabs[0]->show();    DiscSpins[0]->setValue(16); DiscSpins[0]->show();   NDisc++;
    DiscLabs[1]->setText("Axial");          DiscLabs[1]->show();    DiscSpins[1]->setValue(16); DiscSpins[1]->show();  NDisc++;
    DiscLabsExt[0]->setText("Azimuthal");   DiscLabsExt[0]->show(); ExtDiscSpins[0]->setValue(16); ExtDiscSpins[0]->show(); NDiscExt++;
    DiscLabsExt[1]->setText("Axial");       DiscLabsExt[1]->show(); ExtDiscSpins[1]->setValue(16); ExtDiscSpins[1]->show(); NDiscExt++;

    update();
}

void Grid_Options::on_ButtonCylinder_clicked()
{
    // Cylinder type geometry
    GeoType = Half_Cylinder;
    Hide_Labels();
    Hide_Buttons();
    Darken_Button(ui->ButtonCylinder);

    DimLabs[0]->setText("Radius"); DimLabs[0]->show(); DimSpins[0]->setValue(1.0); DimSpins[0]->show(); NDims++;
    DimLabs[1]->setText("Depth");  DimLabs[1]->show(); DimSpins[1]->setValue(10.0); DimSpins[1]->show(); NDims++;
    DimLabsExt[0]->setText("Radial factor"); DimLabsExt[0]->show(); ExtDimSpins[0]->setValue(1.0); ExtDimSpins[0]->show(); NDimsExt++;

    DiscLabs[0]->setText("Azimuthal");      DiscLabs[0]->show();    DiscSpins[0]->setValue(16);  DiscSpins[0]->show();  NDisc++;
    DiscLabs[1]->setText("Radial");         DiscLabs[1]->show();    DiscSpins[1]->setValue(8);   DiscSpins[1]->show();  NDisc++;
    DiscLabs[2]->setText("Vertical");       DiscLabs[2]->show();    DiscSpins[2]->setValue(32);  DiscSpins[2]->show();  NDisc++;
    DiscLabsAux[0]->setText("Radial");      DiscLabsAux[0]->show(); AuxDiscSpins[0]->setValue(8); AuxDiscSpins[0]->show();  NDiscAux++;
    DiscLabsExt[0]->setText("Radial");      DiscLabsExt[0]->show(); ExtDiscSpins[0]->setValue(16); ExtDiscSpins[0]->show();  NDiscExt++;

    update();
}

void Grid_Options::on_ButtonBarge_clicked()
{
    // Barge type geometry

    GeoType = Barge;
    // GeoType = NWT;
    Hide_Labels();
    Hide_Buttons();
    Darken_Button(ui->ButtonBarge);

    DimLabs[NDims]->setText("Length"); DimLabs[NDims]->show(); DimSpins[NDims]->setValue(5.0); DimSpins[NDims]->show(); NDims++;
    DimLabs[NDims]->setText("Width"); DimLabs[NDims]->show(); DimSpins[NDims]->setValue(2.0); DimSpins[NDims]->show(); NDims++;
    DimLabs[NDims]->setText("Draft"); DimLabs[NDims]->show(); DimSpins[NDims]->setValue(2.0); DimSpins[NDims]->show(); NDims++;
    DimLabsExt[NDimsExt]->setText("Length"); DimLabsExt[NDimsExt]->show(); ExtDimSpins[NDimsExt]->setValue(10.0); ExtDimSpins[NDimsExt]->show(); NDimsExt++;
    DimLabsExt[NDimsExt]->setText("Width"); DimLabsExt[NDimsExt]->show(); ExtDimSpins[NDimsExt]->setValue(8.0); ExtDimSpins[NDimsExt]->show(); NDimsExt++;

    DiscLabs[NDisc]->setText("X");          DiscLabs[NDisc]->show();    DiscSpins[NDisc]->setValue(16);  DiscSpins[NDisc]->show();  NDisc++;
    DiscLabs[NDisc]->setText("Y");          DiscLabs[NDisc]->show();    DiscSpins[NDisc]->setValue(8);  DiscSpins[NDisc]->show();  NDisc++;
    DiscLabs[NDisc]->setText("Z");          DiscLabs[NDisc]->show();    DiscSpins[NDisc]->setValue(8);  DiscSpins[NDisc]->show();  NDisc++;
    DiscLabsExt[NDiscExt]->setText("X");    DiscLabsExt[NDiscExt]->show(); ExtDiscSpins[NDiscExt]->setValue(12); ExtDiscSpins[NDiscExt]->show();  NDiscExt++;
    DiscLabsExt[NDiscExt]->setText("Y");    DiscLabsExt[NDiscExt]->show(); ExtDiscSpins[NDiscExt]->setValue(12); ExtDiscSpins[NDiscExt]->show();  NDiscExt++;

    update();
}

void Grid_Options::on_ButtonSparBuoy_clicked()
{
    // Tapered spar buoy type geometry
    GeoType = Tapered_SparBuoy;
    Hide_Labels();
    Hide_Buttons();
    Darken_Button(ui->ButtonSparBuoy);

    DimLabs[NDims]->setText("Lower Radius"); DimLabs[NDims]->show(); DimSpins[NDims]->setValue(6.0); DimSpins[NDims]->show(); NDims++;
    DimLabs[NDims]->setText("Upper radius"); DimLabs[NDims]->show(); DimSpins[NDims]->setValue(4.0); DimSpins[NDims]->show(); NDims++;
    DimLabs[NDims]->setText("Draft"); DimLabs[NDims]->show(); DimSpins[NDims]->setValue(120.0); DimSpins[NDims]->show(); NDims++;
    DimLabs[NDims]->setText("Lower section height"); DimLabs[NDims]->show(); DimSpins[NDims]->setValue(80.0); DimSpins[NDims]->show(); NDims++;
    DimLabs[NDims]->setText("Taper section height"); DimLabs[NDims]->show(); DimSpins[NDims]->setValue(20.0); DimSpins[NDims]->show(); NDims++;
    DimLabsExt[NDimsExt]->setText("Radial factor"); DimLabsExt[NDimsExt]->show(); ExtDimSpins[NDimsExt]->setValue(2.0); ExtDimSpins[NDimsExt]->show(); NDimsExt++;

    DiscLabs[NDisc]->setText("Azimuthal");      DiscLabs[NDisc]->show();    DiscSpins[NDisc]->setValue(16);  DiscSpins[NDisc]->show();  NDisc++;
    DiscLabs[NDisc]->setText("Radial");      DiscLabs[NDisc]->show();    DiscSpins[NDisc]->setValue(8);  DiscSpins[NDisc]->show();  NDisc++;
    DiscLabs[NDisc]->setText("Vertical 1");      DiscLabs[NDisc]->show();    DiscSpins[NDisc]->setValue(40);  DiscSpins[NDisc]->show();  NDisc++;
    DiscLabs[NDisc]->setText("Vertical 2");      DiscLabs[NDisc]->show();    DiscSpins[NDisc]->setValue(10);  DiscSpins[NDisc]->show();  NDisc++;
    DiscLabs[NDisc]->setText("Vertical 3");      DiscLabs[NDisc]->show();    DiscSpins[NDisc]->setValue(10);  DiscSpins[NDisc]->show();  NDisc++;
    DiscLabsAux[NDiscAux]->setText("Radial");      DiscLabsAux[NDiscAux]->show(); AuxDiscSpins[NDiscAux]->setValue(8); AuxDiscSpins[NDiscAux]->show();  NDiscAux++;
    DiscLabsExt[NDiscExt]->setText("Radial");      DiscLabsExt[NDiscExt]->show(); ExtDiscSpins[NDiscExt]->setValue(16); ExtDiscSpins[NDiscExt]->show();  NDiscExt++;

    update();
}

void Grid_Options::on_ButtonOC3_clicked()
{
    // OC3 geometry (dimensions specified)
    GeoType = OC3_SparBuoy;
    Hide_Labels();
    Hide_Buttons();
    Darken_Button(ui->ButtonOC3);

    DimLabs[NDims]->setText("Scaling factor"); DimLabs[NDims]->show(); DimSpins[NDims]->setValue(1.0); DimSpins[NDims]->show(); NDims++;
    DimLabsExt[NDimsExt]->setText("Radial factor"); DimLabsExt[NDimsExt]->show(); ExtDimSpins[NDimsExt]->setValue(1.0); ExtDimSpins[NDimsExt]->show(); NDimsExt++;

    DiscLabs[NDisc]->setText("Azimuthal");      DiscLabs[NDisc]->show();    DiscSpins[NDisc]->setValue(16);  DiscSpins[NDisc]->show();  NDisc++;
    DiscLabs[NDisc]->setText("Radial");      DiscLabs[NDisc]->show();    DiscSpins[NDisc]->setValue(8);  DiscSpins[NDisc]->show();  NDisc++;
    DiscLabs[NDisc]->setText("Vertical 1");      DiscLabs[NDisc]->show();    DiscSpins[NDisc]->setValue(40);  DiscSpins[NDisc]->show();  NDisc++;
    DiscLabs[NDisc]->setText("Vertical 2");      DiscLabs[NDisc]->show();    DiscSpins[NDisc]->setValue(10);  DiscSpins[NDisc]->show();  NDisc++;
    DiscLabs[NDisc]->setText("Vertical 3");      DiscLabs[NDisc]->show();    DiscSpins[NDisc]->setValue(10);  DiscSpins[NDisc]->show();  NDisc++;
    DiscLabsAux[NDiscAux]->setText("Radial");      DiscLabsAux[NDiscAux]->show(); AuxDiscSpins[NDiscAux]->setValue(8); AuxDiscSpins[NDiscAux]->show();  NDiscAux++;
    DiscLabsExt[NDiscExt]->setText("Radial");      DiscLabsExt[NDiscExt]->show(); ExtDiscSpins[NDiscExt]->setValue(16); ExtDiscSpins[NDiscExt]->show();  NDiscExt++;

    update();
}

void Grid_Options::on_ButtonTripleSpar_clicked()
{
    // OC4 geometry (dimensions specified)
    GeoType = Triple_Spar;
    Hide_Labels();
    Hide_Buttons();
    Darken_Button(ui->ButtonTripleSpar);

    DimLabs[NDims]->setText("Radius"); DimLabs[NDims]->show(); DimSpins[NDims]->setValue(5.0); DimSpins[NDims]->show(); NDims++;
    DimLabs[NDims]->setText("Depth"); DimLabs[NDims]->show(); DimSpins[NDims]->setValue(20.0); DimSpins[NDims]->show(); NDims++;
    DimLabs[NDims]->setText("Distance to leg"); DimLabs[NDims]->show(); DimSpins[NDims]->setValue(20.0); DimSpins[NDims]->show(); NDims++;

    DimLabsExt[NDimsExt]->setText("Radial factor"); DimLabsExt[NDimsExt]->show(); ExtDimSpins[NDimsExt]->setValue(1.0); ExtDimSpins[NDimsExt]->show(); NDimsExt++;

    DiscLabs[NDisc]->setText("Azimuthal");      DiscLabs[NDisc]->show();    DiscSpins[NDisc]->setValue(16);  DiscSpins[NDisc]->show();  NDisc++;
    DiscLabs[NDisc]->setText("Radial");         DiscLabs[NDisc]->show();    DiscSpins[NDisc]->setValue(8);  DiscSpins[NDisc]->show();  NDisc++;
    DiscLabs[NDisc]->setText("Vertical");       DiscLabs[NDisc]->show();    DiscSpins[NDisc]->setValue(20);  DiscSpins[NDisc]->show();  NDisc++;

    DiscLabsAux[NDiscAux]->setText("Radial");   DiscLabsAux[NDiscAux]->show(); AuxDiscSpins[NDiscAux]->setValue(8); AuxDiscSpins[NDiscAux]->show();  NDiscAux++;

    DiscLabsExt[NDiscExt]->setText("Azimuthal");    DiscLabsExt[NDiscExt]->show(); ExtDiscSpins[NDiscExt]->setValue(16); ExtDiscSpins[NDiscExt]->show();  NDiscExt++;
    DiscLabsExt[NDiscExt]->setText("Radial");       DiscLabsExt[NDiscExt]->show(); ExtDiscSpins[NDiscExt]->setValue(16); ExtDiscSpins[NDiscExt]->show();  NDiscExt++;

    CheckBoxes[0]->setChecked(false);    CheckBoxes[0]->setText("Turbine \n over \n spar"); CheckBoxes[0]->show();

    update();
}

void Grid_Options::on_ButtonOC4_clicked()
{
    // OC4 geometry (dimensions specified)
    GeoType = SemiSub_OC4;
    Hide_Labels();
    Hide_Buttons();
    Darken_Button(ui->ButtonOC4);

    DimLabs[NDims]->setText("Scaling factor"); DimLabs[NDims]->show(); DimSpins[NDims]->setValue(1.0); DimSpins[NDims]->show(); NDims++;
    DimLabsExt[NDimsExt]->setText("Radial factor"); DimLabsExt[NDimsExt]->show(); ExtDimSpins[NDimsExt]->setValue(1.0); ExtDimSpins[NDimsExt]->show(); NDimsExt++;

    DiscLabs[NDisc]->setText("Azimuthal");      DiscLabs[NDisc]->show();    DiscSpins[NDisc]->setValue(16);  DiscSpins[NDisc]->show();  NDisc++;
    DiscLabs[NDisc]->setText("Radial");         DiscLabs[NDisc]->show();    DiscSpins[NDisc]->setValue(8);  DiscSpins[NDisc]->show();  NDisc++;
    DiscLabs[NDisc]->setText("Vertical 1");      DiscLabs[NDisc]->show();    DiscSpins[NDisc]->setValue(5);  DiscSpins[NDisc]->show();  NDisc++;
    DiscLabs[NDisc]->setText("Vertical 2");      DiscLabs[NDisc]->show();    DiscSpins[NDisc]->setValue(5);  DiscSpins[NDisc]->show();  NDisc++;
    DiscLabs[NDisc]->setText("Vertical 3");      DiscLabs[NDisc]->show();    DiscSpins[NDisc]->setValue(20);  DiscSpins[NDisc]->show();  NDisc++;
    DiscLabsAux[NDiscAux]->setText("Radial");      DiscLabsAux[NDiscAux]->show(); AuxDiscSpins[NDiscAux]->setValue(8); AuxDiscSpins[NDiscAux]->show();  NDiscAux++;
    DiscLabsExt[NDiscExt]->setText("Azimuthal");      DiscLabsExt[NDiscExt]->show(); ExtDiscSpins[NDiscExt]->setValue(16); ExtDiscSpins[NDiscExt]->show();  NDiscExt++;
    DiscLabsExt[NDiscExt]->setText("Radial");      DiscLabsExt[NDiscExt]->show(); ExtDiscSpins[NDiscExt]->setValue(16); ExtDiscSpins[NDiscExt]->show();  NDiscExt++;

    update();
}

void Grid_Options::on_ButtonWigley_clicked()
{
    // Wigley hull geometry
    GeoType = Wigley_Hull;
    Hide_Labels();
    Hide_Buttons();
    Darken_Button(ui->ButtonWigley);

    DimLabs[NDims]->setText("Length");  DimLabs[NDims]->show(); DimSpins[NDims]->setValue(5.0); DimSpins[NDims]->show(); NDims++;
    DimLabs[NDims]->setText("Width");   DimLabs[NDims]->show(); DimSpins[NDims]->setValue(1.0); DimSpins[NDims]->show(); NDims++;
    DimLabs[NDims]->setText("Depth");   DimLabs[NDims]->show(); DimSpins[NDims]->setValue(1.0); DimSpins[NDims]->show(); NDims++;

    DiscLabs[NDisc]->setText("X");      DiscLabs[NDisc]->show();    DiscSpins[NDisc]->setValue(16);  DiscSpins[NDisc]->show();  NDisc++;
    DiscLabs[NDisc]->setText("Z");      DiscLabs[NDisc]->show();    DiscSpins[NDisc]->setValue(16);  DiscSpins[NDisc]->show();  NDisc++;

    update();
}

//--- Generate geometries

void Grid_Options::Generate_Ellipsoid()
{
    // Submerged Ellipsoid is generated
    Boundary = new BEMUse::Ellipsoid();

    Parameters.push_back(BEMUse::Parameter("Semiaxis_a",Real(DimSpins[0]->value())));
    Parameters.push_back(BEMUse::Parameter("Semiaxis_b",Real(DimSpins[1]->value())));
    Parameters.push_back(BEMUse::Parameter("Semiaxis_c",Real(DimSpins[2]->value())));
    Parameters.push_back(BEMUse::Parameter("Ellipsoid_Depth",Real(DimSpins[3]->value())));
    Parameters.push_back(BEMUse::Parameter("NPanels_Axial",DiscSpins[1]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_Azimuthal",DiscSpins[0]->value()));
    if (ui->CosineDist->isChecked())        Parameters.push_back(BEMUse::Parameter("Cosine_Disc",true));
    else                                    Parameters.push_back(BEMUse::Parameter("Cosine_Disc",false));
    // Parameters.push_back(BEMUse::Parameter("Triangular_Panels",true));  // Flag generates triangular panels
}

void Grid_Options::Generate_SemiEllipsoid()
{
    // Submerged Ellipsoid is generated
    Boundary = new BEMUse::Semi_Ellipsoid();

    // Parameters
    Parameters.push_back(BEMUse::Parameter("Semiaxis_a",Real(DimSpins[0]->value())));
    Parameters.push_back(BEMUse::Parameter("Semiaxis_b",Real(DimSpins[1]->value())));
    Parameters.push_back(BEMUse::Parameter("Semiaxis_c",Real(DimSpins[2]->value())));
    Parameters.push_back(BEMUse::Parameter("FreeSurface_Radius",Real(ExtDimSpins[0]->value())));
    Parameters.push_back(BEMUse::Parameter("NPanels_Axial",DiscSpins[1]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_Azimuthal",DiscSpins[0]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_FreeSurface_Radial",ExtDiscSpins[1]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_FreeSurface_Int_Radial",AuxDiscSpins[0]->value()));
    if (ui->CosineDist->isChecked())        Parameters.push_back(BEMUse::Parameter("Cosine_Disc",true));
    else                                    Parameters.push_back(BEMUse::Parameter("Cosine_Disc",false));
    // Parameters.push_back(BEMUse::Parameter("Triangular_Panels",true)); // Flag generates triangular panels
    Parameters.push_back(BEMUse::Parameter("Triangular_Panels",true)); // Flag generates triangular panels
}

void Grid_Options::Generate_Cylinder()
{
    Boundary = new BEMUse::Half_Cylinder();

    // Parameters
    Parameters.push_back(BEMUse::Parameter("Radius", Real(DimSpins[0]->value())));
    Parameters.push_back(BEMUse::Parameter("Draft", Real(DimSpins[1]->value())));
    Parameters.push_back(BEMUse::Parameter("FreeSurface_Radius",Real(ExtDimSpins[0]->value())));
    Parameters.push_back(BEMUse::Parameter("NPanels_Azimuthal",DiscSpins[0]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_Radial",DiscSpins[1]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_Vertical",DiscSpins[2]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_FreeSurface_Radial",ExtDiscSpins[0]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_FreeSurface_Int_Radial",AuxDiscSpins[0]->value()));
    if (ui->CosineDist->isChecked())        Parameters.push_back(BEMUse::Parameter("Cosine_Disc",true));
    else                                    Parameters.push_back(BEMUse::Parameter("Cosine_Disc",false));
}

void Grid_Options::Generate_Barge()
{
    Boundary = new BEMUse::Barge();

    // Parameters
    Parameters.push_back(BEMUse::Parameter("Tank_Length", Real(DimSpins[0]->value())));
    Parameters.push_back(BEMUse::Parameter("Tank_Width", Real(DimSpins[1]->value())));
    Parameters.push_back(BEMUse::Parameter("Tank_Depth",Real(DimSpins[2]->value())));
    Parameters.push_back(BEMUse::Parameter("NPanels_Tank_Length",DiscSpins[0]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_Tank_Width",DiscSpins[1]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_Tank_Depth",DiscSpins[2]->value()));
    Parameters.push_back(BEMUse::Parameter("FreeSurface_Length",Real(ExtDimSpins[0]->value())));
    Parameters.push_back(BEMUse::Parameter("FreeSurface_Width",Real(ExtDimSpins[1]->value())));
    Parameters.push_back(BEMUse::Parameter("NPanels_FreeSurface_Length",ExtDiscSpins[0]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_FreeSurface_Width",ExtDiscSpins[1]->value()));
    if (ui->CosineDist->isChecked())        Parameters.push_back(BEMUse::Parameter("Cosine_Disc",true));
    else                                    Parameters.push_back(BEMUse::Parameter("Cosine_Disc",false));
}

void Grid_Options::Generate_SparBuoy()
{
    Boundary = new BEMUse::Tapered_SparBuoy();

    // Parameters
    Parameters.push_back(BEMUse::Parameter("Lower_Radius", Real(DimSpins[0]->value())));
    Parameters.push_back(BEMUse::Parameter("Upper_Radius", Real(DimSpins[1]->value())));
    Parameters.push_back(BEMUse::Parameter("Draft", Real(DimSpins[2]->value())));
    Parameters.push_back(BEMUse::Parameter("Height_Taper_Begin", Real(DimSpins[3]->value())));
    Parameters.push_back(BEMUse::Parameter("Height_Taper_End", Real(DimSpins[4]->value())));
    Parameters.push_back(BEMUse::Parameter("FreeSurface_Radius", Real(ExtDimSpins[0]->value())));
    Parameters.push_back(BEMUse::Parameter("NPanels_Radial", DiscSpins[1]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_Vertical_Section1", DiscSpins[2]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_Vertical_Section2", DiscSpins[3]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_Vertical_Section3", DiscSpins[4]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_Azimuthal", DiscSpins[0]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_FreeSurface_Radial", ExtDiscSpins[0]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_FreeSurface_Int_Radial", AuxDiscSpins[0]->value()));
    if (ui->CosineDist->isChecked())        Parameters.push_back(BEMUse::Parameter("Cosine_Disc",true));
    else                                    Parameters.push_back(BEMUse::Parameter("Cosine_Disc",false));
}

void Grid_Options::Generate_TripleSpar()
{
    Boundary = new BEMUse::Triple_Spar();

    // Parameters
    Parameters.push_back(BEMUse::Parameter("Radius", Real(DimSpins[0]->value())));
    Parameters.push_back(BEMUse::Parameter("Draft", Real(DimSpins[1]->value())));
    Parameters.push_back(BEMUse::Parameter("Spar_Radius", Real(DimSpins[2]->value())));
    Parameters.push_back(BEMUse::Parameter("FreeSurface_Radius",Real(ExtDimSpins[0]->value())));
    Parameters.push_back(BEMUse::Parameter("NPanels_Azimuthal",DiscSpins[0]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_Radial",DiscSpins[1]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_Vertical",DiscSpins[2]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_FreeSurface_Radial",ExtDiscSpins[0]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_FreeSurface_Int_Radial",AuxDiscSpins[0]->value()));
    if (ui->CosineDist->isChecked())        Parameters.push_back(BEMUse::Parameter("Cosine_Disc",true));
    else                                    Parameters.push_back(BEMUse::Parameter("Cosine_Disc",false));
    Parameters.push_back(BEMUse::Parameter("Shift_Over_Leg",CheckBoxes[0]->isChecked()));
}

void Grid_Options::Generate_OC3SparBuoy()
{
    Boundary = new BEMUse::OC3_SparBuoy();

    // Parameters
    Parameters.push_back(BEMUse::Parameter("Scaling_Factor", Real(DimSpins[0]->value())));
    Parameters.push_back(BEMUse::Parameter("FreeSurface_Radius", Real(ExtDimSpins[0]->value())));
    Parameters.push_back(BEMUse::Parameter("NPanels_Radial", DiscSpins[1]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_Vertical_Section1", DiscSpins[2]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_Vertical_Section2", DiscSpins[3]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_Vertical_Section3", DiscSpins[4]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_Azimuthal", DiscSpins[0]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_FreeSurface_Radial", ExtDiscSpins[0]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_FreeSurface_Int_Radial", AuxDiscSpins[0]->value()));
    if (ui->CosineDist->isChecked())        Parameters.push_back(BEMUse::Parameter("Cosine_Disc",true));
    else                                    Parameters.push_back(BEMUse::Parameter("Cosine_Disc",false));
    Parameters.push_back(BEMUse::Parameter("Triangular_Panels", false));
}

void Grid_Options::Generate_OC4TripleSpar()
{
    Boundary = new BEMUse::SemiSub_OC4();

    // Parameters

    // Spar Leg
    Parameters.push_back(BEMUse::Parameter("Scaling_Factor",Real(DimSpins[0]->value())));
    Parameters.push_back(BEMUse::Parameter("NPanels_Azimuthal",DiscSpins[0]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_Radial_Base",DiscSpins[1]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_Vertical_HeavePlate",DiscSpins[2]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_Horizontal_HeavePlate",DiscSpins[3]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_Vertical",DiscSpins[4]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_FreeSurface_Int_Radial",AuxDiscSpins[0]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_FreeSurface_Radial",DiscSpins[1]->value()));
    Parameters.push_back(BEMUse::Parameter("FreeSurface_Radius",Real(ExtDimSpins[0]->value())));

    // Central column
    Parameters.push_back(BEMUse::Parameter("NPanels_Radial",DiscSpins[1]->value()));

    if (ui->CosineDist->isChecked())        Parameters.push_back(BEMUse::Parameter("Cosine_Disc",true));
    else                                    Parameters.push_back(BEMUse::Parameter("Cosine_Disc",false));
    Parameters.push_back(BEMUse::Parameter("Triangular_Panels",false)); // Generate quadratic elements for this config
    Parameters.push_back(BEMUse::Parameter("Shift_Over_Leg",CheckBoxes[0]->isChecked())); // Generate quadratic elements for this config
}

void Grid_Options::Generate_WigleyHull()
{
    Boundary = new BEMUse::Wigley_Hull();

    Parameters.push_back(BEMUse::Parameter("Length", Real(DimSpins[0]->value())));
    Parameters.push_back(BEMUse::Parameter("Beam",  Real(DimSpins[1]->value())));
    Parameters.push_back(BEMUse::Parameter("Draft", Real(DimSpins[2]->value())));
    Parameters.push_back(BEMUse::Parameter("NPanels_Length",DiscSpins[0]->value()));
    Parameters.push_back(BEMUse::Parameter("NPanels_Depth",DiscSpins[1]->value()));
    if (ui->CosineDist->isChecked())        Parameters.push_back(BEMUse::Parameter("Cosine_Disc",true));
    else                                    Parameters.push_back(BEMUse::Parameter("Cosine_Disc",false));
}

