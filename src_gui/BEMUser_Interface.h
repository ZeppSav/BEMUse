#ifndef BEMUSER_INTERFACE_H
#define BEMUSER_INTERFACE_H

//--------------------
//--- Qt Includes ----
//--------------------

#include <QOpenGLWidget>
#include <QThread>
#include <QMouseEvent>
#include <QApplication>
#include <QTimer>
#include <QElapsedTimer>

//--------------------
//---- Dialogues -----
//--------------------

#include "Grid_Options.h"
#include "Solver_Setup.h"

//--------------------
//--- BEMUse Source---
//--------------------

#include "Boundary/Boundary_Base.h"
#include "Solver/Solver_Base.h"
#include "Visualise.h"

#define GLPANELS                1

enum Status_Type
{
    kSuccess,           // No problems
    kH_Mesh_Above_Zero  // There are mesh elements above the free surface line
};

static std::vector<std::string> Status_Message =
{
    "No errors detected",                                                               // kSuccess
    "Mesh elements have been detected above the free surface. Check your mesh!"         // kH_Mesh_Above_Zero
};

extern Status_Type Status;

namespace Ui {class BEMUser_Interface;}

// Creating a QObject class to handle executing solver in a separate thread.

static Real WheelRollFactor = 120;
static Real zoomFactor = 1.1;
static Real TimeFactor = TwoPI/3000;

class WorkerThread : public QThread
{
    Q_OBJECT

public:

    //--- Solver object
    BEMUse::Solver *Solver = NULL;
    BEMUse::Boundary *Boundary = NULL;

    //--- Simulation arrays
    StdVector FreqList;
    StdVector BetaList;

    void run() override
    {
        QString result;
        timestamp_t t0 = get_timestamp();               // Begin clock
        for (int i=0; i<FreqList.size(); i++)
        {
            Solver->Set_Real(FreqList[i]);
            Solver->Solve();
            analysis_progress(i+1,FreqList.size());
        }
        emit analysis_complete(result);

        // Store visualisation objects
        StdAppend(Boundary->RadSolArray,Solver->RadSolArray);
        StdAppend(Boundary->DiffSolArray,Solver->DiffSolArray);
        StdAppend(Boundary->FS_Rad_Array,Solver->FS_Rad_SolArray);
        StdAppend(Boundary->FS_Scat_Array,Solver->FS_Scat_SolArray);

        Real Tsec = (get_timestamp() - t0) / 1000.0L;   // Record end time
        std::cout << "Execution complete. Time taken: " << Tsec << " ms." << std::endl;

        // Call finished in order to delete this thread.
        emit finished();
    }

//    ~WorkerThread()  {}     // This is to avoid issues with deleting the object

signals:
    void analysis_progress(int n, int ntot);
    void analysis_complete(const QString &s);
    void finished();
};

class BEMUser_Interface : public QOpenGLWidget
{
    Q_OBJECT

private slots:

    void on_Button_Grid_clicked();                  // Generate grid
    void on_Configure_clicked();                    // Solver configuration
    void on_Button_Execute_clicked();               // Execute solver

private:

    Ui::BEMUser_Interface *ui;

    //--- BEMUse architecture
    BEMUse::Boundary *Boundary = NULL;
    BEMUse::Solver *Solver = NULL;

    //--- Simulation arrays
    StdVector FreqList;
    StdVector BetaList;

    //--- Solver Thread
    WorkerThread *workerThread = NULL;

    //--- Visualisation Vars
    bool Sol_Avail = false;
    int FView = 0;      // Which frequency do we want to visualise
    int FTot = 0;       // What is the maximum number of frequencies?
    int DOFView = 0;    // Which DOF do we want to view?
    int BetaView = 0;    // Which DOF do we want to view?
    int BetaTot = 0;    // Which DOF do we want to view?

    bool ShowNormals = false;
    bool ShowPhiJ = false;
    bool ShowPhiD = false;
    bool ShowSurfacepanels = true;
    bool ShowInteriorFreeSurface = false;
    bool ShowExteriorFreeSurface = false;

    // Free surface visulisation vars
    bool ShowFreeSurface = false;
    bool ShowScatteredFreeSurface = false;

    QTimer *Vis_Trigger;
    QElapsedTimer *Vis_Timer;

    //---- Position and scaling vars
    int xRot;                               // X rotation view angle
    int yRot;                               // Y rotation view angle
    int zRot;                               // Z rotation view angle
    QPoint lastPos;                         // Relative mouse position
    float HeightScale;
    float WidthScale;
    float zoomScale = 1.0;                  // Zoomscaling
    float Axes_Scale = 1.0;                 // Coordinate axes scale
    GLfloat xshift=0, yshift=0, zshift=0;   // Translation variables

    //---- Drawing functions
    void draw_coordinate_axes();
    void draw_scene();
    void Clear_Vis();

private slots:      // Visualisation actions

    void Update_Status(std::string M = "");         // Status updates (text prompt)

    //--- Buttons
    void on_NormVis_toggled(bool checked)           {ShowNormals = checked;}
    void on_BoundaryVis_toggled(bool checked)       {ShowSurfacepanels = checked;}
    void on_IFSVis_toggled(bool checked)            {ShowInteriorFreeSurface = checked;}
    void on_EFSVis_toggled(bool checked)            {ShowExteriorFreeSurface = checked;}
    void on_ExportButton_clicked();
    void on_ImportButton_clicked();
    void on_GeoDelete_clicked();
    void on_SurfPot_toggled(bool checked);
    void on_SurfDiff_toggled(bool checked);
    void on_WaveVis_toggled(bool checked);
    void on_WaveVis_2_toggled(bool checked);

public:
    BEMUser_Interface(QOpenGLWidget *parent = 0);
    ~BEMUser_Interface();

public slots:
    void setXRotation(int angle);
    void setYRotation(int angle);
    void setZRotation(int angle);
    void updateResults(const int &i1, const int &i2)    {Update_Status("Analysis: " + std::to_string(i1) + " of " + std::to_string(i2) + " complete.");}
    void handleResults(const QString &)                 {Update_Status("Hydrodynamic analysis complete."); Sol_Avail = true;}

    //--- Slots for external modules
    void Receive_Geometry(BEMUse::Boundary *B);
    void Receive_Solver(BEMUse::Solver *S, StdVector Freqs, StdVector Betas);

signals:
    void xRotationChanged(int angle);
    void yRotationChanged(int angle);
    void zRotationChanged(int angle);
    void geometryChanged();
    void operate(const QString &);

protected:

    //--- Virtual QOpenGLWidget functions
    void initializeGL();
    void paintGL();
    void resizeGL(int width, int height);
    void mouseDoubleClickEvent(QMouseEvent *event);
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void wheelEvent(QWheelEvent *event);

};

//inline void QAppend(QString &S, QString A)                 {S.append(A);}
//inline void QFillAppend(QString &S, char *A, int NMax)     {QAppend(S,QString("%1").arg(A, -NMax, QChar(' ')));}
//inline void QFillAppend(QString &S, QString A, int NMax)   {QAppend(S,QString("%1").arg(A, -NMax, QChar(' ')));}
//inline void QFillAppend(QString &S, Real A, int NMax)      {QFillAppend(S,QString::number(A,'f',2),NMax);}
//inline void QFillAppendLong(QString &S, Real A, int NMax)  {QFillAppend(S,QString("%1").arg(A,4,'e',5),NMax);}

#endif // BEMUSER_INTERFACE_H
