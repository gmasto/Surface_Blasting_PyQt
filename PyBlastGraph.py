"""
PyBlastGraph is a program that produces 2 or 3 dimensional representations of a blast site
"""
from PyQt5 import QtCore, QtWidgets, uic
from pyqtgraph import PlotWidget
import pyqtgraph as pg
import sys
import numpy as np
import pyqtgraph.opengl as gl
from PyQt5.QtGui import QVector3D
import BlastGraph_Inputs as Dp


class MainWindow(QtWidgets.QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        # Load the UI Page
        uic.loadUi('PyBlastGraph.ui', self)
        # Set the window title
        self.setWindowTitle("PyBlastGraph")

        # Connect the radiobuttons, checkboxes and buttons with their corresponding methods
        self.radioButton_3D.toggled.connect(self.on_dimension_selection)
        self.radioButton_2D.toggled.connect(self.on_dimension_selection)

        self.radioButton_1Face.toggled.connect(self.on_geometry_selection)
        self.radioButton_2Face.toggled.connect(self.on_geometry_selection)

        self.presplitting_checkBox.stateChanged.connect(self.on_presplit_selection)
        self.stag_checkBox.stateChanged.connect(self.on_pattern_selection)

        self.design_pushButton.clicked.connect(self.on_design_clicked)

        # Set the value of the program's components when it executes for the first time
        self.radioButton_3D.setChecked(True)
        self.dimensions = 3

        self.radioButton_1Face.setChecked(True)
        self.geo2 = "n"

        self.presplitting_checkBox.setChecked(False)
        self.presplit = "n"

        self.stag_checkBox.setChecked(False)
        self.stag = "n"

        # the widget for the 3D graph
        self.openGLWidget.setVisible(True)
        # the widget for the 2D graph
        self.plotwidget.setVisible(False)
        gl.GLViewWidget.setBackgroundColor(self.openGLWidget, 'k')

    # The methods below institute the functionality of the program
    def on_dimension_selection(self):
        """
        The on_dimension_selection method gets the selection of the user from the Graph Dimensions radio buttons group
        and controls the visibility of the graph widgets

        :return:
        """
        value = self.sender().text()
        if value == "3-D":
            self.dimensions = 3
            self.openGLWidget.setVisible(True)
            self.plotwidget.setVisible(False)
        else:
            self.dimensions = 2
            self.openGLWidget.setVisible(False)
            self.plotwidget.setVisible(True)

    def on_geometry_selection(self):
        """
        The on_geometry_selection method gets the selection of the user from the Free Face Selection radio buttons group
        and sets the program to calculate the results accordingly.
        :return:
        """
        value = self.sender().text()
        if value == "1 Free Face":
            self.geo2 = "n"
        else:
            self.geo2 = "y"

    def on_presplit_selection(self, state: bool):
        """
        The on_presplit_selection method checks the state of the "Design Presplit Holes" checkbox and sets the
        program to calculate the results accordingly.
        :param state: bool, the state of the checkbox
        :return: None
        """
        if state == QtCore.Qt.Checked:
            self.presplit = "y"
        else:
            self.presplit = "n"

    def on_pattern_selection(self, state: bool):
        """
        The on_pattern_selection method checks the state of the "Apply Staggered Blast Pattern" checkbox and sets the
        program to calculate the results accordingly.
        :param state: bool, the state of the checkbox
        :return: None
        """
        if state == QtCore.Qt.Checked:
            self.stag = "y"
        else:
            self.stag = "n"

    def on_design_clicked(self):
        """
        The on_design_clicked method executes when the user clicks the "Design" button. Performs the calculation needed
        by the program to display the graphs in the graph widgets according to the previous selections by the user.
        :return:
        """
        print("Design Button Clicked")
        # Clear old plots before adding the new ones
        for item in reversed(self.openGLWidget.items):
            self.openGLWidget.removeItem(item)

        self.plotwidget.clear()

        # Initiations for suppressing weak warnings (not really needed)
        xpresplit = None
        ypresplit = None
        xpresplitb = None
        ypresplitb = None
        hps = None
        hops = None
        dhopsc = None
        dhpsc = None

        # Insert rows and holes per row
        rows = int(Dp.rows)
        perrow = int(Dp.perrow)

        # Insert parameters
        b = Dp.bp                                   # the burden
        s = Dp.sp                                   # the spacing
        a = 360 - (90 - Dp.a)                       # the face angle
        u = Dp.u                                    # the subdrill
        h = Dp.h                                    # the height of the blast hole
        hb = Dp.hb                                  # the length of the bottom charge
        ho = Dp.ho                                  # the length of the stemming
        hc = Dp.hc                                  # the length of the column charge
        dz = np.cos(2 * np.pi - a * np.pi / 180)    # the horizontal correction
        dyc = np.sin(2 * np.pi - a * np.pi / 180)   # the vertical correction
        hcc = hc / dz                               # the blast's hole's column height correction
        hbc = hb / dz                               # the blast's hole's bottom height correction
        dhcc = (hc - hcc)                           # additional corrections to the height
        dhbc = (hb - hbc)                           # additional corrections to the height

        # check the free face selection
        if self.geo2 == "y":
            a1 = Dp.a1  # the angle between the free faces of the blast site
        else:
            a1 = 90     # if 1 free face is selected

        # Coordinates corrections along the X axis
        dxb = b * np.tan(np.pi / 2 - a1 * np.pi / 180)
        dxs2 = (s / 2) * np.tan(np.pi / 2 - a1 * np.pi / 180)

        # Create the staggered coordinates of the holes (X, Y)
        if self.stag == "y":
            perrow += 1
            x = np.zeros((rows, perrow))
            y = np.zeros((rows, perrow))
            for j in range(0, rows):
                if j % 2 != 0:
                    de = np.zeros(perrow)
                    for i in range(2, perrow):
                        if self.geo2 == "n":
                            de[0] = s / 2
                        else:
                            de[0] = b

                        de[1] = de[0] + s / 2
                        de[i] = de[i - 1] + s

                    de[-1] = de[-2] + s / 2
                    x[j, :] = de
                else:
                    do = np.zeros(perrow)
                    for i in range(len(do)):
                        if self.geo2 == "n":
                            do[i] = (s / 2) + i * s
                        else:
                            do[i] = b + i * s

                    do[-1] = 0
                    x[j, :] = do

            if self.geo2 == "y":
                for i in range(rows):
                    x[i, :] = x[i, :] - (i * dxb)
                    if i % 2 == 0:
                        x[i, -1] = 0

            for j in range(rows):
                if j % 2 != 0:
                    re = np.zeros(perrow)
                    for i in range(len(re)):
                        re[i] = (j + 1) * b

                    y[j, :] = re
                else:
                    ro = np.zeros(perrow)
                    for i in range(len(ro)):
                        ro[i] = (j + 1) * b

                    ro[-1] = 0
                    y[j, :] = ro

        else:
            dholes = np.zeros(perrow)
            for i in range(1, perrow):
                if self.geo2 == "y":
                    dholes[0] = b - dxb
                else:
                    dholes[0] = s / 2

                dholes[i] = dholes[i - 1] + s

            rholes = np.zeros(rows)
            for i in range(rows):
                rholes[i] = (i + 1) * b

            x, y = np.meshgrid(dholes, rholes)

            if a1 != 90:
                for i in range(1, rows):
                    x[i, :] = x[i - 1, :] - dxb

        if self.stag == "y":
            perrow = perrow - 1

        dlimits = np.zeros(perrow + 2)
        for i in range(2, perrow + 1):
            if self.geo2 == "y":
                dlimits[1] = b
            else:
                dlimits[1] = s / 2

            dlimits[i] = dlimits[i - 1] + s
            dlimits[-1] = dlimits[-2] + s / 2

        extralims = 2

        if self.geo2 == "y":
            rlimits = np.zeros(rows + 2 + extralims)
            for i in range(1, rows + 2 + extralims):
                rlimits[i] = i * b

            rlimits[-1] = rlimits[-2] + s / 2

        else:
            rlimits = np.zeros(rows + 1 + extralims)
            for i in range(len(rlimits)):
                rlimits[i] = i * b

        xlimits, ylimits = np.meshgrid(dlimits, rlimits)

        if a1 != 90:
            for i in range(1, rows + 1 + extralims):
                xlimits[i, :] = xlimits[i - 1, :] - dxb

            xlimits[-1, :] = xlimits[-2, :] - dxs2

        if self.stag == "y":
            xlimits += dxb

        # PRESPLITTING
        if self.presplit == "y":
            hops = Dp.hops
            hps = Dp.hps
            hopsc = hops / dz
            hpsc = hps / dz
            dhopsc = (hops - hopsc)
            dhpsc = (hps - hpsc)
            sps = Dp.sps
            dysps = sps - sps * np.cos(np.pi / 2 - a1 * np.pi / 180)
            rsps = sps * np.cos(np.pi / 2 - a1 * np.pi / 180)
            xsps = rsps * np.tan(np.pi / 2 - a1 * np.pi / 180)

            if self.geo2 == "n":
                xpresplit = np.ones([int(Dp.nps), 2])
                xpresplit[:, 0] = xpresplit[:, 0] * xlimits[0, 0]
                xpresplit[:, 1] = xpresplit[:, 1] * xlimits[0, -1]
                ypresplit = np.ones([int(Dp.nps), 2])
                dpr = np.zeros([1, int(Dp.nps)]).T
                for i in range(1, int(Dp.nps) + 1):
                    dpr[i - 1] = rows * Dp.bp - (i - 1) * sps

                ypresplit = ypresplit * dpr
            else:
                # Side presplitting Holes
                xpresplit = np.ones([int(Dp.nps2), 1])
                xpresplit = xpresplit * xlimits[-1-extralims, -1]
                for i in range(1, int(Dp.nps2)):
                    xpresplit[i] = xpresplit[i - 1] + xsps

                ypresplit = np.ones([int(Dp.nps2), 1])
                dpr = np.zeros([1, int(Dp.nps2)]).T
                for i in range(1, int(Dp.nps2) + 1):
                    dpr[i - 1] = rows * b + (s / 2) - (i - 1) * (sps - dysps)

                ypresplit = ypresplit * dpr
                # Back presplitting Holes
                ypresplitb = np.ones([1, int(Dp.npsback2)]) * ylimits[-1-extralims, -1]
                ypresplitb = ypresplitb.T
                xpresplitb = np.ones([1, int(Dp.npsback2)]).T * xlimits[-1-extralims, -1] - sps
                for i in range(1, int(Dp.npsback2)):
                    xpresplitb[i] = xpresplitb[i - 1] - sps

        xminlim = x.min(initial=None)
        xmaxlim = x.max(initial=None)
        yminlim = y.min(initial=None)
        ymaxlim = y.max(initial=None)
        xlimmax = xlimits.max()
        xlimmin = xlimits.min()

        if self.stag == "y":
            perrow = perrow + 1

        # 3D Holes and Surfaces
        if self.dimensions == 3:
            self.openGLWidget.setCameraPosition(pos=QVector3D(xmaxlim/2, ymaxlim/2, 0),
                                                distance=2*xlimmax, azimuth=225)
            # cylinder (Blast Holes)
            for j in range(rows):
                for i in range(perrow):
                    if (y[j, i] and x[j, i]) != 0:
                        md = gl.MeshData.cylinder(rows=10, cols=20, radius=[0.08, 0.08], length=hb)
                        m5 = gl.GLMeshItem(meshdata=md, color=(1, 0, 0, 1), smooth=False, drawEdges=False,
                                           edgeColor=(1, 0, 0, 1), shader='balloon')
                        m5.rotate(a, a, 0, 0)
                        m5.translate(x[j, i], y[j, i], 0)
                        self.openGLWidget.addItem(m5)

                        md = gl.MeshData.cylinder(rows=10, cols=20, radius=[0.08, 0.08], length=hc)
                        m6 = gl.GLMeshItem(meshdata=md, color=(0, 0, 1, 1), smooth=False, drawEdges=False,
                                           edgeColor=(1, 0, 0, 1), shader='balloon')
                        m6.rotate(a, a, 0, 0)
                        m6.translate(x[j, i], y[j, i] + (hb * dyc), hb+dhbc)
                        self.openGLWidget.addItem(m6)

                        md = gl.MeshData.cylinder(rows=10, cols=20, radius=[0.08, 0.08], length=ho)
                        m7 = gl.GLMeshItem(meshdata=md, color=(0.22, 1, 0, 1), smooth=False, drawEdges=False,
                                           edgeColor=(1, 0, 0, 1), shader='balloon')
                        m7.rotate(a, a, 0, 0)
                        m7.translate(x[j, i], y[j, i] + (hb * dyc + hc * dyc), (hb+hc)+(dhbc+dhcc))
                        self.openGLWidget.addItem(m7)

            # Presplit Holes
            if self.presplit == "y":
                rx, cx = xpresplit.shape
                xpresplit = xpresplit.reshape(rx * cx)
                ypresplit = ypresplit.reshape(rx * cx)
                if self.geo2 == "n":
                    for j in range(int(rx * cx)):
                        md = gl.MeshData.cylinder(rows=10, cols=20, radius=[0.035, 0.035], length=hps)
                        m5 = gl.GLMeshItem(meshdata=md, color=(1, 1, 0, 1), smooth=False, drawEdges=False,
                                           edgeColor=(1, 0, 0, 1), shader='balloon')
                        m5.rotate(a, a, 0, 0)
                        m5.translate(xpresplit[j], ypresplit[j], u)
                        self.openGLWidget.addItem(m5)

                        md = gl.MeshData.cylinder(rows=10, cols=20, radius=[0.035, 0.035], length=hops)
                        m6 = gl.GLMeshItem(meshdata=md, color=(0.22, 1, 0, 1), smooth=False, drawEdges=False,
                                           edgeColor=(1, 0, 0, 1), shader='balloon')
                        m6.rotate(a, a, 0, 0)
                        m6.translate(xpresplit[j], ypresplit[j] + (hps * dyc), u-dhopsc+hps+dhpsc)
                        self.openGLWidget.addItem(m6)
                else:
                    rxb, cxb = xpresplitb.shape
                    xpresplitb = xpresplitb.reshape(rxb * cxb)
                    ypresplitb = ypresplitb.reshape(rxb * cxb)
                    # Side PS holes
                    for j in range(int(rx * cx)):
                        md = gl.MeshData.cylinder(rows=10, cols=20, radius=[0.035, 0.035], length=hps)
                        m5 = gl.GLMeshItem(meshdata=md, color=(1, 1, 0, 1), smooth=False, drawEdges=False,
                                           edgeColor=(1, 0, 0, 1), shader='balloon')
                        m5.rotate(a, a, 0, 0)
                        m5.translate(xpresplit[j], ypresplit[j], u)
                        self.openGLWidget.addItem(m5)

                        md = gl.MeshData.cylinder(rows=10, cols=20, radius=[0.035, 0.035], length=hops)
                        m6 = gl.GLMeshItem(meshdata=md, color=(0.22, 1, 0, 1), smooth=False, drawEdges=False,
                                           edgeColor=(1, 0, 0, 1), shader='balloon')
                        m6.rotate(a, a, 0, 0)
                        m6.translate(xpresplit[j], ypresplit[j, ] + (hps * dyc), u-dhopsc+hps+dhpsc)
                        self.openGLWidget.addItem(m6)

                    # Back PS holes
                    for j in range(int(rxb * cxb)):
                        md = gl.MeshData.cylinder(rows=10, cols=20, radius=[0.035, 0.035], length=hps)
                        m5 = gl.GLMeshItem(meshdata=md, color=(1, 1, 0, 1), smooth=False, drawEdges=False,
                                           edgeColor=(1, 0, 0, 1), shader='balloon')
                        m5.rotate(a, a, 0, 0)
                        m5.translate(xpresplitb[j], ypresplitb[j], u)
                        self.openGLWidget.addItem(m5)

                        md = gl.MeshData.cylinder(rows=10, cols=20, radius=[0.035, 0.035], length=hops)
                        m6 = gl.GLMeshItem(meshdata=md, color=(0.22, 1, 0, 1), smooth=False, drawEdges=False,
                                           edgeColor=(1, 0, 0, 1), shader='balloon')
                        m6.rotate(a, a, 0, 0)
                        m6.translate(xpresplitb[j], ypresplitb[j] + (hps * dyc), u-dhopsc+hps+dhpsc)
                        self.openGLWidget.addItem(m6)

            if self.geo2 == "y":
                x1 = xlimits[0, 0]
                x2 = xlimits[-1, 0]
                x3 = x1
                x6 = x2
            else:
                x1 = xlimmin-10
                x2 = xlimmin-10
                x3 = x1
                x6 = x2

            verts = np.array([                                      # #            6 _________ 5
                [xlimmax+10, 0, u],                                 # 0             |\         \
                [x1, 0, u],                                         # 1             | \         \
                [x2, ylimits[-1, 0]+h * dyc, u],                    # 2           2 |  \_________\ 7
                [x3, h * dyc, h * dz],                              # 3             \  |3        |<---------4
                [xlimmax+10, ylimits[-1, 0]+h * dyc, u],            # 4              \ |         |
                [xlimmax+10, ylimits[-1, 0]+h * dyc, h * dz],       # 5               \|_________|
                [x6, ylimits[-1, 0]+h * dyc, h * dz],               # 6                1         0
                [xlimmax+10, h * dyc, h * dz]                       # 7
            ])
            # Create triangular faces (2 for each 'square' face). The numbers are the peaks of the edges.
            if self.geo2 == "y":
                faces = np.array([
                    [1, 0, 7], [1, 3, 7],
                    [1, 2, 4], [1, 0, 4],
                    [1, 2, 6], [1, 3, 6],
                    [3, 6, 5], [3, 7, 5]
                ])
            else:
                faces = np.array([
                    [1, 0, 7], [1, 3, 7],
                    [1, 2, 4], [1, 0, 4],
                    [3, 6, 5], [3, 7, 5]
                ])

            colors = np.array([[0.79, 0.46, 0.17, 0.3] for i in range(8)])  # in range(No of faces)
            sf = gl.GLMeshItem(vertexes=verts, faces=faces, faceColors=colors, smooth=False)
            sf.setGLOptions('additive')
            self.openGLWidget.addItem(sf)

            # Vertices and face for the bottom surface
            verts = np.array([                              # #            2___________3
                [xlimmax + 10, -10, u],                     # 0             \          \
                [xlimmin-10, -10, u],                       # 1              \          \
                [xlimmin-10, ylimits[-1, 0]+h * dyc, u],    # 2               \__________\
                [xlimmax + 10, ylimits[-1, 0]+h * dyc, u],  # 3               1           0
            ])
            faces = np.array([
                [1, 0, 3], [1, 2, 3]
            ])
            colors = np.array([[0.79, 0.46, 0.17, 0.3] for i in range(2)])  # in range(No of faces)
            sb = gl.GLMeshItem(vertexes=verts, faces=faces, faceColors=colors, smooth=False)
            sb.setGLOptions('additive')
            self.openGLWidget.addItem(sb)

            # Limits Line Plot
            if self.geo2 == "y":
                pv = -3
            else:
                pv = -1
            # Top line limits
            pts = np.array([
                [xlimits[0, 0], h * dyc, h * dz],
                [xlimits[pv, 0], ylimits[-3, 0]+h * dyc, h * dz],
                [xlimits[pv, -1], ylimits[-3, 0]+h * dyc,  h * dz],
                [xlimits[0, -1], h * dyc, h * dz],
            ])
            linelims = gl.GLLinePlotItem(pos=pts, color=(1, 0, 1, 1), width=0.1, antialias=False)
            self.openGLWidget.addItem(linelims)

            # Back line limits
            pts = np.array([
                [xlimits[pv, 0], ylimits[-3, 0]+h * dyc, h * dz],
                [xlimits[pv, 0], ylimits[-3, 0], u],
                [xlimits[pv, -1], ylimits[-3, 0]+h * dyc,  h * dz],
                [xlimits[pv, -1], ylimits[-3, 0],  u],
            ])
            linelims = gl.GLLinePlotItem(pos=pts, color=(1, 0, 1, 1), width=0.1, antialias=False)
            self.openGLWidget.addItem(linelims)

            # Down line limits
            pts = np.array([
                [xlimits[0, 0], ylimits[0, 0], u],
                [xlimits[pv, 0], ylimits[-3, 0], u],
                [xlimits[pv, -1], ylimits[-3, 0],  u],
                [xlimits[0, -1], ylimits[0, 0], u],
            ])
            linelims = gl.GLLinePlotItem(pos=pts, color=(1, 0, 1, 1), width=0.1, antialias=False)
            self.openGLWidget.addItem(linelims)

        # 2D Holes and Limits Plots
        else:
            PlotWidget.setBackground(self.plotwidget, 'w')
            pen = pg.mkPen(color=(255, 0, 0), style=QtCore.Qt.NoPen)
            legbrush = pg.mkBrush(color=(150, 150, 150), style=QtCore.Qt.NoPen)
            legend = pg.LegendItem(offset=(0., .5), brush=legbrush, labelTextColor=(0, 0, 0))
            legend.setParentItem(self.plotwidget.getPlotItem())
            x = x.reshape(rows*perrow)
            y = y.reshape(rows*perrow)
            x = np.delete(x, np.where(y == 0))
            y = np.delete(y, np.where(y == 0))
            pl1 = self.plotwidget.plot(x, y,
                                       name="holes", pen=pen, symbol='o', symbolSize=7, symbolBrush='r')
            legend.addItem(pl1, "holes")

            self.plotwidget.setTitle("Plot Title", color="k", size="30pt")
            styles = {'color': 'k', 'font-size': '20px'}
            self.plotwidget.setLabel('left', "Side", **styles)
            self.plotwidget.setLabel('bottom', "Front", **styles)
            self.plotwidget.setXRange(xminlim - 10, xmaxlim + 10, padding=0)
            self.plotwidget.setYRange(yminlim - 10, ymaxlim + 10, padding=0)

            if self.presplit == "y":
                pen = pg.mkPen(color=(0, 0, 0), style=QtCore.Qt.NoPen)
                if self.geo2 == "n":
                    for hole in range(len(xpresplit)):
                        self.plotwidget.plot(xpresplit[hole, :], ypresplit[hole, :],
                                             pen=pen, symbol='o', symbolSize=4, symbolBrush='k')
                        if hole == (len(xpresplit)-1):
                            pl2 = self.plotwidget.plot(xpresplit[hole, :], ypresplit[hole, :], name="presplit holes",
                                                       pen=pen, symbol='o', symbolSize=4, symbolBrush='k')
                            legend.addItem(pl2, "presplit holes")

                else:
                    for hole in range(len(xpresplit)):
                        self.plotwidget.plot(xpresplit[hole, :], ypresplit[hole, :],
                                             pen=pen, symbol='o', symbolSize=4, symbolBrush='k')
                    for hole in range(len(xpresplitb)):
                        self.plotwidget.plot(xpresplitb[hole, :], ypresplitb[hole, :],
                                             pen=pen, symbol='o', symbolSize=4, symbolBrush='k')
                        if hole == (len(xpresplitb)-1):
                            pl3 = self.plotwidget.plot(xpresplitb[hole, :], ypresplitb[hole, :], name="presplit holes",
                                                       pen=pen, symbol='o', symbolSize=4, symbolBrush='k')
                            legend.addItem(pl3, "presplit holes")

            pen = pg.mkPen(color=(0, 0, 255), style=QtCore.Qt.SolidLine)

            xbox = [xlimits[0, 0], xlimits[-3, 0], xlimits[-3, -1], xlimits[0, -1], xlimits[0, 0]]
            ybox = [ylimits[0, 0], ylimits[-3, 0], ylimits[-3, -1], ylimits[0, -1], ylimits[0, 0]]
            pl4 = self.plotwidget.plot(xbox, ybox, name="limits", pen=pen)

            legend.addItem(pl4, "limits")

        print("Complete")

        # # Prints
        # print(f"a = {a}")
        # print(f"dy = {dyc}")
        # print(f"dx = {dz}")
        # print(f"tan(a1) = {dxl}")
        # print("=" * 100)
        # print(f"x coordinates = \n{x}")
        # print("=" * 100)
        # print(f"y coordinates = \n{y}")
        # print("=" * 100)
        # if self.presplit == "y":
        #     print(xpresplit)
        #     print(ypresplit)
        #     print(dpr)
        #     print(f"xpscoor = {xpresplit}")
        #     print(f"ypscoor = {ypresplit}")
        #     print(f"ypscoorB = {ypresplitB}")
        #     print(f"xpscoorB = {xpresplitB}")
        # print("=" * 100)
        # print(f"x limits = \n{xlimits}")
        # print("=" * 100)
        # print(f"y limits = \n{ylimits}")
        # print("=" * 100)
        # print(f"min = {xminlim}")
        # print(f"max = {xmaxlim}")
        # print("=" * 100)
        # print(f"limplot = {xlimits[0, -1]+10}")
        # print("=" * 100)
        # print(f"maxlim = {xlimmax}")
        # print(f"minlim = {xlimmin}")
        # print(f"dxs2 = {dxs2}")
        # print("=" * 100)
        # if self.dimensions == 2:
        #     print(xbox)
        #     print(ybox)


def main():
    app = QtWidgets.QApplication(sys.argv)
    mainwindow = MainWindow()
    mainwindow.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
