#include "BEMUser_Interface.h"
#include "ui_bemuser_interface.h"

#include <QDebug>
#include <QApplication>
#include <QFileDialog>
#include <GL/glu.h>

Status_Type Status;
std::string Output_Status = "";

#include <sstream>
template <typename T>
std::string StdStringPrec(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

BEMUser_Interface::BEMUser_Interface(QOpenGLWidget  *parent) :  QOpenGLWidget(parent), ui(new Ui::BEMUser_Interface)
{
    ui->setupUi(this);
}

//--- Virtual QOpenGLWidget functions ----

void BEMUser_Interface::initializeGL()
{
    // Set background color
    glClearColor(1,1,1,1);      // Set white background

    // Set initial oriential angle
    xRot = -45*16;
    yRot = 0;
    zRot = 45*16;

    // Initialize trigger & connect slot
    Vis_Trigger = new QTimer(this);
    connect(Vis_Trigger, SIGNAL(timeout()), this, SLOT(update()));
    Vis_Trigger->start(100);

    // Initialize timer
    Vis_Timer = new QElapsedTimer();
    Vis_Timer->start();
}

void BEMUser_Interface::paintGL()
{
    // Clear previous lists

//    glEnable(GL_BLEND);                                     // Disable for opaque elements
//    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);      // Disable for opaque elements
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_POINT_SMOOTH);
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    glDisable(GL_MULTISAMPLE);
    glHint(GL_MULTISAMPLE_FILTER_HINT_NV, GL_NICEST);
    glEnable(GL_DEPTH_TEST);

    glLoadIdentity();
    glScaled(zoomScale, zoomScale, zoomScale);
    glTranslatef(xshift, yshift, -10.0);
    glRotatef(xRot / 16.0, 1.0, 0.0, 0.0);
    glRotatef(yRot / 16.0, 0.0, 1.0, 0.0);
    glRotatef(zRot / 16.0, 0.0, 0.0, 1.0);

    draw_scene();               // Generate scene
    draw_coordinate_axes();     // Coordinate system axes

    // Display lists
    glCallList(GL_GREATER);

    // Update output box
//    ui->Output_Label->setText(QString::fromStdString(Output_Status));
}

void BEMUser_Interface::resizeGL(int width, int height)
{
    int side = qMin(width, height);
    int range = 10;
    //glViewport((width - side) / 2, (height - side) / 2, side, side);
    glViewport(0, 0, width, height);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    glOrtho(-0.5*range*zoomScale, +0.5*range*zoomScale, -0.5*height/width*range*zoomScale, +0.5*height/width*range*zoomScale, -5*range, +5*range);

    glMatrixMode(GL_MODELVIEW);
}

//--- Drawing functions ----

void BEMUser_Interface::draw_coordinate_axes()
{
    // This function visualises the coordinate axes of the scene

    glColor4d(0,0,0,1);
    glLineWidth(2);

    // X axis
    glBegin(GL_LINES);
    glVertex3f(0,0,0);
    glVertex3f(Axes_Scale,0,0);
    glEnd();

    // Y axis
    glBegin(GL_LINES);
    glVertex3f(0,0,0);
    glVertex3f(0,Axes_Scale,0);
    glEnd();


    // Z axis
    glBegin(GL_LINES);
    glVertex3f(0,0,0);
    glVertex3f(0,0,Axes_Scale);
    glEnd();
}

void BEMUser_Interface::draw_scene()
{
    //---

    if (Boundary==NULL) return;

//    if (glIsList(GLPANELS))  glDeleteLists(GLPANELS, 1); // Panel

    //--- Create drawing elements

    if (ShowNormals)                    BEMUse::Visualise_Centroid_CS(Boundary);
    if (ShowSurfacepanels)              BEMUse::Visualise_Boundary_Elements(Boundary);
    if (ShowInteriorFreeSurface)        BEMUse::Visualise_Interior_FreeSurface_Elements(Boundary);
    if (ShowExteriorFreeSurface)        BEMUse::Visualise_Exterior_FreeSurface_Elements(Boundary);
    if (Sol_Avail && ShowPhiJ)          BEMUse::Visualise_Radiation_Solution(Boundary, FView, DOFView);
    if (Sol_Avail && ShowPhiD)          BEMUse::Visualise_Diffraction_Solution(Boundary, FView, BetaView);

    glEnable(GL_BLEND);                                     // Disable for opaque elements
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);      // Disable for opaque elements
    if (Sol_Avail && ShowFreeSurface)           BEMUse::Visualise_Free_Surface(Boundary, FView, DOFView, Real(Vis_Timer->elapsed())*TimeFactor);
    if (Sol_Avail && ShowScatteredFreeSurface)  BEMUse::Visualise_ScatteredFree_Surface(Boundary, FView, BetaView, Real(Vis_Timer->elapsed())*TimeFactor);

    glEndList();

    //--- Update output
    ui->Output_Label->setText(QString::fromStdString(Output_Status));
}

//--- Mouse events ----

void BEMUser_Interface::mousePressEvent(QMouseEvent *event)
{
    lastPos = event->pos();
}

void BEMUser_Interface::wheelEvent(QWheelEvent *event)
{
    // Note: A single "roll" of the wheel corresponds to approx delta y of 120.
    // We should therefore scale appropriately...

//    int DeltaF = int (-1.0*event->angleDelta().y()/WheelRollFactor);

    if (QGuiApplication::queryKeyboardModifiers().testFlag(Qt::ControlModifier) && Sol_Avail)
    {
        //--- Toggle visualisation frequency
        // Note: A single "roll" of the wheel corresponds to approx delta y of 120.
        // We should therefore scale appropriately...
        int DeltaF = int (-1.0*event->angleDelta().y()/WheelRollFactor);
        int FNew = FView + DeltaF;
        if      (FNew>=FTot)    FView = FTot-1;
        else if (FNew<0)        FView = 0;
        else                    FView = FNew;

        // Update message
        Update_Status("Visualising solution DOF " + std::to_string(DOFView+1) + " at frequency " + StdStringPrec(FreqList[FView],3) + " Hz.");
        update();
    }
    else if (QGuiApplication::queryKeyboardModifiers().testFlag(Qt::ShiftModifier) && Sol_Avail)
    {
        //--- Toggle visualisation degree of freedom
        int DeltaF = int (-1.0*event->angleDelta().y()/WheelRollFactor);
        int DOFNew = DOFView + DeltaF;
        if      (DOFNew>=5)     DOFView = 5;
        else if (DOFNew<0)      DOFView = 0;
        else                    DOFView = DOFNew;

        // Update message
        Update_Status("Visualising solution DOF " + std::to_string(DOFView+1) + " at frequency " + StdStringPrec(FreqList[FView],3) + "Hz.");
        update();
    }
    else if (QGuiApplication::queryKeyboardModifiers().testFlag(Qt::AltModifier) && Sol_Avail)
    {
        //--- Toggle visualisation wave angle
        int DeltaF = int (-1.0*event->angleDelta().x()/WheelRollFactor);
        int BetaNew = BetaView + DeltaF;
        if      (BetaNew>=BetaTot)  BetaView = BetaTot-1;
        else if (BetaNew<0)         BetaView = 0;
        else                        BetaView = BetaNew;

        // Update message
        Update_Status("Visualising solution DOF " + std::to_string(DOFView+1) + " at wave angle " + StdStringPrec(BetaList[BetaView]*R2D,2) + " degrees.");
        update();
    }
    else
    {
        //--- Zoom in & out
        if (event->angleDelta().y() < 0)  zoomScale /= zoomFactor;
        if (event->angleDelta().y() > 0)  zoomScale *= zoomFactor;

        GLint viewport[4];
        glGetIntegerv(GL_VIEWPORT,viewport);
        resizeGL(viewport[2], viewport[3]);

        update(); // call paintGL()
    }
}

void BEMUser_Interface::mouseMoveEvent(QMouseEvent *event)
{
    if(event->buttons() == Qt::RightButton)
    {
        // Translate
        int dx = event->x() - lastPos.x();
        int dy = event->y() - lastPos.y();
        xshift = 0.005/zoomScale*dx;
        yshift = -0.005/zoomScale*dy;
    }
    else
    {
        // Rotate
        int dx = event->x() - lastPos.x();
        int dy = event->y() - lastPos.y();
        setXRotation(xRot + 8 * dy);
        setYRotation(yRot + 8 * dx);
        lastPos = event->pos();
    }

    update();
}

void BEMUser_Interface::mouseDoubleClickEvent(QMouseEvent *event)
{
    // Look down z axis: (0,0,0)
    // Look down y axis: (-90,0,0)
    // Look down z axis: (-90,0,90)
    setXRotation(-45*16);
    setYRotation(0);
    setZRotation(45*16);
}

static void qNormalizeAngle(int &angle)
{
    while (angle < 0)
        angle += 360 * 16;
    while (angle > 360 * 16)
        angle -= 360 * 16;
}

void BEMUser_Interface::setXRotation(int angle)
{
    qNormalizeAngle(angle);
    if (angle != xRot) {
        xRot = angle;
        emit xRotationChanged(angle);
        update();
    }
}

void BEMUser_Interface::setYRotation(int angle)
{
    qNormalizeAngle(angle);
    if (angle != yRot) {
        yRot = angle;
        emit yRotationChanged(angle);
        update();
    }
}

void BEMUser_Interface::setZRotation(int angle)
{
    qNormalizeAngle(angle);
    if (angle != zRot) {
        zRot = angle;
        emit zRotationChanged(angle);
        update();
    }
}

BEMUser_Interface::~BEMUser_Interface()
{
    delete ui;
}

//--- Geometry configuration

void BEMUser_Interface::on_Button_Grid_clicked()
{
    // This generates a geometry for analysis within BEMUse

    if (Boundary!=NULL)
    {
        Update_Status("A geometry has already been specified. Don't be greedy!");
        return;
    }

    // Generate new grid options window
    Grid_Options *GOD = new Grid_Options(0);
    connect(GOD, &Grid_Options::Pass_Geometry, this, &BEMUser_Interface::Receive_Geometry);
    GOD->show();

    Status = kSuccess;

    if (Status!=kSuccess)   {Update_Status(Status_Message[Status]);}
    else                    {Update_Status("If the mesh looks ok, then go ahead and carry out preprocessing.");}

    update();                                                           //--- Repaint scene
}

#define WAMIT_FILES "WAMIT geometry file (*.gdf)"
#define NEMOH_FILES "NEMOH geometry file (*.mar)"
#define STL_FILES   "STL file (*.stl)"
#define PNL_FILES   "PNL file (*.pnl)"

void BEMUser_Interface::on_ImportButton_clicked()
{
    if (Boundary!=NULL) return;
    QString filter = WAMIT_FILES ";;" NEMOH_FILES ";;" STL_FILES ";;" PNL_FILES;
    QString selectedFilter;
    QString Pfad = QFileDialog::getOpenFileName(this, tr("Find File"), QDir::currentPath(),filter,&selectedFilter);
    if (Pfad.isEmpty()) {Update_Status("No File selected.");   return;}
    if (selectedFilter==WAMIT_FILES)    Boundary = new BEMUse::GDF_Geometry();
    if (selectedFilter==NEMOH_FILES)    Boundary = new BEMUse::MAR_Geometry();
    if (selectedFilter==STL_FILES)      Boundary = new BEMUse::STL_Geometry();
    if (selectedFilter==STL_FILES)      Boundary = new BEMUse::PNL_Geometry();

//    Boundary->Read_Input_File(Pfad);
    std::string Pfadstring = Pfad.toStdString();
    Boundary->Read_Input_File(Pfadstring);
    Boundary->Setup();

    if (Status!=kSuccess)   {Update_Status(Status_Message[Status]);}
    else                    {Update_Status("If the mesh looks ok, then go ahead and carry out preprocessing.");}
    if (Boundary!=NULL)     ShowSurfacepanels = true;                   //--- Activate visualisation flag
    update();                                                           //--- Repaint scene
}

void BEMUser_Interface::on_ExportButton_clicked()
{
    if (Boundary==NULL) return;
    QString filter = WAMIT_FILES ";;" NEMOH_FILES;
    QString selectedFilter;
    QString Pfad = QFileDialog::getSaveFileName(this, tr("Save File"), QDir::currentPath(),filter,&selectedFilter);
    std::string StdPfad = Pfad.toStdString();
    if (selectedFilter==WAMIT_FILES)    Boundary->Export_Geometry_GDF(StdPfad);
    if (selectedFilter==NEMOH_FILES)    Boundary->Export_Geometry_MAR(StdPfad);

//    // Debug out external node positions
//    std::vector<BEMUse::SP_Node> ExtNodes;
//    Boundary->Get_Ext_Nodes(ExtNodes);
//    for (BEMUse::SP_Node N : ExtNodes)
//    {
//        Vector3 P = N->Position_Local();
//        std::cout << P(0) << " " << P(1) << " " << P(2) << std::endl;
//     }
}

void BEMUser_Interface::Update_Status(std::string M)
{
    // This is a helper function to update the status and refresh the screen.
    if (!M.empty())    Output_Status = M;
    update();
}

void BEMUser_Interface::on_GeoDelete_clicked()
{
    if (Boundary!=NULL){
        delete Boundary;
        Boundary = NULL;
        Clear_Vis();
    }
    if (Solver!=NULL){
        delete Solver;
        Solver = NULL;
        FreqList.clear();
        BetaList.clear();
        Clear_Vis();
    }
    if (workerThread!=NULL) workerThread = NULL;
}

void BEMUser_Interface::Clear_Vis()
{
    // Clears the visualisation options
    ui->BoundaryVis->setChecked(false);
    ui->IFSVis->setChecked(false);
    ui->EFSVis->setChecked(false);
    ui->NormVis->setChecked(false);
    ui->SurfPot->setChecked(false);
    ui->SurfDiff->setChecked(false);
    ui->WaveVis->setChecked(false);
    ui->WaveVis_2->setChecked(false);

    Sol_Avail = false;
    ShowNormals = false;
    ShowPhiJ = false;
    ShowPhiD = false;
    ShowSurfacepanels = false;
    ShowInteriorFreeSurface = false;
    ShowExteriorFreeSurface = false;

    ShowFreeSurface = false;
    ShowScatteredFreeSurface = false;
}

//--- Solver configuration

void BEMUser_Interface::on_Configure_clicked()
{
    // This prepares the solver

    //--- Error check
    if (Status!=kSuccess)       {Update_Status(Status_Message[Status]);             return;}
    if (Boundary==NULL)         {Update_Status("Geometry must be specified.");      return;}
    if (Solver!=NULL)           {Update_Status("Solver already defined.");          return;}

    // Generate new solver setup dialogue
    Solver_Setup *Setup = new Solver_Setup(0);
    connect(Setup, &Solver_Setup::Pass_Solver, this, &BEMUser_Interface::Receive_Solver);
    Setup->show();
}

void BEMUser_Interface::Receive_Geometry(BEMUse::Boundary *B)
{
    // The geometry created within the geometry creation module is passed.
    Boundary = B;

    // Set vis flags.
    ui->BoundaryVis->setChecked(true);
    ShowSurfacepanels = true;
}

void BEMUser_Interface::Receive_Solver(BEMUse::Solver *S, StdVector Freqs, StdVector Betas)
{
    // The solver setup dialogue allows the user to specify solver parameters.
    // The solver is generated there and is passed it to the main window.
    // Here the additional tasks (such as passing the solver to the external thread to not block up
    // main window rendering) is carried out.
    Solver = S;

    //--- Create new Worker thread object and connect slots
    workerThread = new WorkerThread();
    connect(workerThread, &WorkerThread::analysis_complete, this, &BEMUser_Interface::handleResults);
    connect(workerThread, &WorkerThread::analysis_progress, this, &BEMUser_Interface::updateResults);
    connect(workerThread, &WorkerThread::finished, workerThread,  &QObject::deleteLater);
    connect(workerThread, SIGNAL(finished()), workerThread, SLOT(deleteLater()), Qt::DirectConnection);     // Deletes object upon completion

    //--- Store values there & carry out setup of geometry
    workerThread->Solver = Solver;
    workerThread->Boundary = Boundary;
    workerThread->Solver->Setup(Boundary);

    //--- Solver setup
    StdAppend(FreqList,Freqs);
    StdAppend(BetaList,Betas);
    StdAppend(workerThread->FreqList,FreqList);
    StdAppend(workerThread->BetaList,BetaList);

//    for (int i=0; i<FreqList.size(); i++) std::cout << FreqList[i] << std::endl;
//    for (int i=0; i<BetaList.size(); i++) std::cout << BetaList[i] << std::endl;

    //--- Update output message
    Update_Status("Ok. Preliminaries have worked. Go ahead and carry out analysis.");

    //--- Repaint scene
    update();

    //--- Prepare visualisation vars
    FTot = FreqList.size();
    BetaTot = BetaList.size();

    std::cout << "Solver Configured succesfully: " << FTot << " frequencies, " << BetaTot << " analysis angles." << std::endl;
}

//--- Execute solver

void BEMUser_Interface::on_Button_Execute_clicked()
{
    // Initialise worked thread.
    workerThread->start();
}

//--- Button events

void BEMUser_Interface::on_SurfPot_toggled(bool checked)
{
    if (ShowPhiD){
        ShowPhiD = false;
        ui->SurfDiff->setChecked(false);
    }
    ShowPhiJ = checked;
    ui->SurfPot->setChecked(checked);
}

void BEMUser_Interface::on_SurfDiff_toggled(bool checked)
{
    if (ShowPhiJ){
        ShowPhiJ = false;
        ui->SurfPot->setChecked(false);
    }
    ShowPhiD = checked;
    ui->SurfDiff->setChecked(checked);
}

void BEMUser_Interface::on_WaveVis_toggled(bool checked)
{
    if (ShowScatteredFreeSurface){
        ShowScatteredFreeSurface = false;
        ui->WaveVis_2->setChecked(false);
    }
    ShowFreeSurface = checked;
    ui->WaveVis->setChecked(checked);
}

void BEMUser_Interface::on_WaveVis_2_toggled(bool checked)
{
    if (ShowFreeSurface){
        ShowFreeSurface = false;
        ui->WaveVis->setChecked(false);
    }
    ShowScatteredFreeSurface = checked;
    ui->WaveVis_2->setChecked(checked);
}

