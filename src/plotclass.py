from abc import ABC, abstractmethod
import plotly.graph_objects as go
import numpy as np


class Plot(ABC):
    """Abstract class for create graphs"""

    def check_data_and_layout(self):
        prefix = f"\nclass: '{self.__class__.__name__}', data: '{self._data is not None}', " \
                 f"layout: '{self._layout is not None}'"
        if self._data is None and self._layout is None:
            error_text = prefix + "\nObject has no defined 'data' and 'layout'."
        elif self._data is None:
            error_text = prefix + "\nObject has no defined 'data'."
        else:
            error_text = prefix + "\nObject has no defined 'layout'."

        return error_text

    def object2plot(self):
        """Method for pack object`s data and layout to go.Figure()"""

        if self._data is not None and self._layout is not None:
            return go.Figure(data=self._data, layout=self._layout)
        else:
            raise ValueError(self.check_data_and_layout())


class ModelScatter3D(Plot):
    """Class PlotScatter3D"""

    # default_marker_size = 4
    # default_marker_color = "#85c5ed"  # "#219dad"
    # active_marker_size = 10
    # active_marker_color = '#ed4242'  # "#49d4e6"

    def _set_default_layout(self):
        """Method for set default layout"""

        self._layout = dict(margin={'l': 0, 'r': 0, 't': 0, 'b': 0}, plot_bgcolor='#111111',
                            paper_bgcolor='rgba(0,0,0,0)',
                            hoverdistance=-1,
                            scene=dict(
                                xaxis=dict(nticks=8, title='x', backgroundcolor='rgba(0, 0, 0, 0)', gridcolor='white',
                                           gridwidth=1, color='white'),
                                yaxis=dict(nticks=8, title='y', backgroundcolor='rgba(0, 0, 0, 0)',
                                           gridcolor='white', gridwidth=1, color='white'),
                                zaxis=dict(nticks=8, title='z', backgroundcolor='rgba(0, 0, 0, 0)', gridcolor='white',
                                           gridwidth=1, color='white'),
                                camera=dict(up={'x': 0, 'y': 0, 'z': 1}, center={'x': 0, 'y': 0, 'z': -0.1},
                                            eye={'x': -1, 'y': -1, 'z': 1}),
                                aspectratio={"x": 0.7, "y": 0.7, "z": 0.6}))

    def update_layout(self, new_layout):
        """Set new layout for the PlotScatter3D object"""
        self._layout = new_layout

    def _set_data(self):
        """Set data for the PlotScatter3D object"""
        import pandas as pd

        def calculate_solutions(ab, bc, cd, eq_zero_var):
            stp = 0.1
            min_val = 0
            max_val = 8
            AB = float(ab)
            BC = float(bc)
            CD = float(cd)

            if eq_zero_var in ['mx', 'my', 'nx', 'ny']:

                if eq_zero_var == "mx":
                    nx = ny = np.arange(min_val, max_val, stp)
                    mx, my = 0, np.zeros((len(nx), len(ny)))
                    for i in range(0, len(nx)):
                        for j in range(0, len(ny)):
                            my[i, j] = (-AB * (mx + nx[i] + ny[j]) + 2 * mx * nx[i]) / (-AB - 2 * ny[j])
                    x, y, z = nx, ny, my
                    x_label, y_label, z_label = "nx", "ny", "my"

                elif eq_zero_var == "my":
                    nx = ny = np.arange(min_val, max_val, stp)
                    my, mx = 0, np.zeros((len(nx), len(ny)))
                    for i in range(0, len(nx)):
                        for j in range(0, len(ny)):
                            mx[i, j] = (-AB * (nx[i] + ny[j] - my) + 2 * my * ny[j]) / (AB - 2 * nx[i])
                    x, y, z = nx, ny, mx
                    x_label, y_label, z_label = "nx", "ny", "mx"

                elif eq_zero_var == "nx":
                    mx = my = np.arange(min_val, max_val, stp)
                    nx, ny = 0, np.zeros((len(mx), len(my)))
                    for i in range(0, len(mx)):
                        for j in range(0, len(my)):
                            ny[i, j] = (-AB * (mx[i] + nx - my[j]) + 2 * mx[i] * nx) / (AB - 2 * my[j])
                    x, y, z = mx, my, ny
                    x_label, y_label, z_label = "mx", "my", "ny"

                elif eq_zero_var == "ny":
                    mx = my = np.arange(min_val, max_val, stp)
                    ny, nx = 0, np.zeros((len(mx), len(my)))
                    for i in range(0, len(mx)):
                        for j in range(0, len(my)):
                            nx[i, j] = (-AB * (mx[i] + ny - my[j]) + 2 * my[j] * ny) / (AB - 2 * mx[i])
                    x, y, z = mx, my, nx
                    x_label, y_label, z_label = "mx", "my", "nx"

            elif eq_zero_var in ['px', 'py', 'sx', 'sy']:

                if eq_zero_var == "px":
                    sx = sy = np.around(np.arange(min_val, max_val, stp), 1)
                    px, py = 0, np.zeros((len(sx), len(sy)))
                    for i in range(0, len(sx)):
                        for j in range(0, len(sy)):
                            py[i, j] = (BC * (px - sx[i] - sy[j]) - 2 * sx[i] * px) / (-BC - 2 * sy[j])
                    x, y, z = sx, sy, py
                    x_label, y_label, z_label = "sx", "sy", "py"

                elif eq_zero_var == "py":
                    sx = sy = np.around(np.arange(min_val, max_val, stp), 1)
                    py, px = 0, np.zeros((len(sx), len(sy)))
                    for i in range(0, len(sx)):
                        for j in range(0, len(sy)):
                            px[i, j] = (BC * (-sx[i] - sy[j] + py) + 2 * sy[j] * py) / (-BC + 2 * sx[i])
                    x, y, z = sx, sy, px
                    x_label, y_label, z_label = "sx", "sy", "px"

                elif eq_zero_var == "sx":
                    px = py = np.around(np.arange(min_val, max_val, stp), 1)
                    sx, sy = 0, np.zeros((len(px), len(py)))
                    for i in range(0, len(px)):
                        for j in range(0, len(py)):
                            sy[i, j] = (BC * (px[i] - sx + py[j]) - 2 * sx * px[i]) / (BC - 2 * py[j])
                    x, y, z = px, py, sy
                    x_label, y_label, z_label = "px", "py", "sy"

                elif eq_zero_var == "sy":
                    px = py = np.around(np.arange(min_val, max_val, stp), 1)
                    sy, sx = 0, np.zeros((len(px), len(py)))
                    for i in range(0, len(px)):
                        for j in range(0, len(py)):
                            sx[i, j] = (BC * (px[i] - sy + py[j]) + 2 * sy * py[j]) / (BC + 2 * px[i])
                    x, y, z = px, py, sx
                    x_label, y_label, z_label = "px", "py", "sx"

            return x, x_label, y, y_label, z, z_label

        def find_decimals(value):
            import decimal
            return abs(decimal.Decimal(str(value)).as_tuple().exponent)

        x, self.__x_l, y, self.__y_l, z, self.__z_l = calculate_solutions(self.ab, self.bc, self.cd, self.eq_zero_var)
        self.__plot = pd.DataFrame(columns=['x', 'y', 'z'])
        for i in range(0, np.size(x)):
            for j in range(0, np.size(y)):
                if find_decimals(z[i, j]) == 1:
                    new_row = {'x': x[i], 'y': y[j], 'z': z[i, j]}
                    self.__plot = self.__plot.append(new_row, ignore_index=True)

        tags = []
        plot_as_np = self.__plot.to_numpy()
        xl = self.__x_l[0] + '<sub>' + self.__x_l[1] + '</sub>'
        yl = self.__y_l[0] + '<sub>' + self.__y_l[1] + '</sub>'
        zl = self.__z_l[0] + '<sub>' + self.__z_l[1] + '</sub>'

        for row in plot_as_np:
            tag = '<b><i>' + xl + ':</b></i> ' + ' ' * (4 - len("{0:.1f}".format(row[0]))) + \
                  "{0:.1f}".format(row[0]) + ' [mm]<br>' + \
                  '<b><i>' + yl + ':</b></i> ' + ' ' * (4 - len("{0:.1f}".format(row[1]))) + \
                  "{0:.1f}".format(row[1]) + ' [mm]<br>' + \
                  '<b><i>' + zl + ':</b></i> ' + ' ' * (4 - len("{0:.1f}".format(row[2]))) + \
                  "{0:.1f}".format(row[2]) + ' [mm]' + '<extra></extra>'
            tags.append(tag)

        self.__prepared_data = dict(
            x=np.array(self.__plot['x']),
            y=np.array(self.__plot['y']),
            z=np.array(self.__plot['z']),
            type='scatter3d',
            hoverlabel=dict(align='left', bgcolor='#103d52', bordercolor='#ffd64f', font=dict(family="Consolas")),
            hovertemplate = tags,
            mode="markers",
            marker=dict(
                opacity=1,
                size=[self._default_marker_size for k in range(len(self.__plot['x']))],
                color=[self._default_marker_color for k in range(len(self.__plot['x']))],
                line=dict(width=8, color='#1573ad'))
        )
        self._data = go.Scatter3d(self.__prepared_data)

    def _input_validation(self):
        """Method for validation inputs"""

        conditionals_list = [
            self.__ab_range[0] <= self.ab <= self.__ab_range[1],
            self.__bc_range[0] <= self.bc <= self.__bc_range[1],
            self.__cd_range[0] <= self.cd <= self.__cd_range[1],
            self.eq_zero_var in self.__eq_zero_var_range]

        warnings_list = [
            f"\n[ab]: value '{self.ab}' is out of range. 'ab' must be in "
            f"range <{self.__ab_range[0]}-{self.__ab_range[1]}>",
            f"\n[bc]: value '{self.bc}' is out of range. 'bc' must be in "
            f"range <{self.__bc_range[0]}-{self.__bc_range[1]}>",
            f"\n[cd]: value '{self.cd}' is out of range. 'cd' must be in "
            f"range <{self.__cd_range[0]}-{self.__cd_range[1]}>",
            f"\n[eq_zero_var]: name '{self.eq_zero_var}' is out of range. 'eq_zero_var' must be from "
            f"the following list: {self.__eq_zero_var_range}"]

        warnings = ""

        for idx, conditional in enumerate(conditionals_list):
            if not conditional:
                warnings += warnings_list[idx]
                # TODO destructor when any conditional in conditionals_list returned False
        if len(warnings) > 0:
            raise ValueError(warnings)

    def __init__(self):
        self._default_marker_size = 4
        self._default_marker_color = "#85c5ed"
        self._hovered_marker_size = 10
        self._hovered_marker_color = "#1dcfb1"
        self._clicked_marker_size = 12
        self._clicked_marker_color = "#67cf1d"
        self._rejected_marker_size = 8
        self._rejected_marker_color = "#a32c2c"

    def build(self, ab: int, bc: int, cd: int, eq_zero_var: str):
        # self._default_marker_color = PlotScatter3D.default_marker_color
        # self._default_marker_size = PlotScatter3D.default_marker_size
        self.__ab_range = [24, 50]
        self.__bc_range = [24, 50]
        self.__cd_range = [12, 37]
        self.__eq_zero_var_range = ['mx', 'my', 'nx', 'ny', 'px', 'py', 'sx', 'sy']
        self.ab = ab
        self.bc = bc
        self.cd = cd
        self.eq_zero_var = eq_zero_var
        self._input_validation()
        self._set_data()
        self._set_default_layout()
        self.full_figure = self.object2plot()

    def update_hovered_markers(self, full_figure, hd, cd1):

        markers = full_figure["data"][0]["marker"]
        markers["size"] = [self._default_marker_size for k in markers["size"]]
        markers["color"] = [self._default_marker_color for k in markers["color"]]

        if cd1 is not None:
            markers["size"][cd1['points'][0]['pointNumber']] = self._clicked_marker_size
            markers["color"][cd1['points'][0]['pointNumber']] = self._clicked_marker_color
        if hd is not None:
            markers["size"][hd['points'][0]['pointNumber']] = self._hovered_marker_size
            markers["color"][hd['points'][0]['pointNumber']] = self._hovered_marker_color

        full_figure["data"][0]["marker"] = markers
        return full_figure

    def update_clicked_markers(self, full_figure, cd):

        x = full_figure["data"][0]["x"]
        y = full_figure["data"][0]["y"]
        z = full_figure["data"][0]["z"]

        x_clicked = cd["points"][0]["x"]
        y_clicked = cd["points"][0]["y"]
        z_clicked = cd["points"][0]["z"]

        markers = full_figure["data"][0]["marker"]
        markers["size"] = [self._default_marker_size for k in markers["size"]]
        markers["color"] = [self._default_marker_color for k in markers["color"]]

        for i in range(len(x)):
            if (x_clicked == x[i]) or (y_clicked == y[i]) or (z_clicked == z[i]):
                markers["size"][i] = self._rejected_marker_size
                markers["color"][i] = self._rejected_marker_color


        full_figure["data"][0]["marker"] = markers
        return full_figure



class MyFigure(go.Figure):

    default_marker_size = 4
    default_marker_color = "#85c5ed"
    clicked_marker_size = 12
    clicked_marker_color = "#67cf1d"
    hovered_marker_size = 10
    hovered_marker_color = "#1dcfb1"

    def update_hovered_markers(self, hd, cd1):

        markers = self.data[0]["marker"]
        markers["size"] = [self.default_marker_size for k in markers["size"]]
        markers["color"] = [self.default_marker_color for k in markers["color"]]
        markers_size = list(markers["size"])
        markers_color = list(markers["color"])

        if cd1 is not None:
            markers_size[cd1['points'][0]['pointNumber']] = self.clicked_marker_size
            markers_color[cd1['points'][0]['pointNumber']] = self.clicked_marker_color

        if hd is not None:
            markers_size[hd['points'][0]['pointNumber']] = self.hovered_marker_size
            markers_color[hd['points'][0]['pointNumber']] = self.hovered_marker_color

        markers["size"] = tuple(markers_size)
        markers["color"] = tuple(markers_color)
        self.data[0]["marker"] = markers


class FigureDict:

    def __init__(self, data, layout):
        # self._data = fdict["data"]
        # self._layout = fdict["layout"]
        # self._fdict = fdict
        self.data = data
        self.layout = layout

    def update_hovered_markers(self, hd, marker_styles):

        default_marker_size = 4
        default_marker_color = "#85c5ed"
        clicked_marker_size = 12
        clicked_marker_color = "#67cf1d"
        hovered_marker_size = 10
        hovered_marker_color = "#1dcfb1"

        markers = self.data[0]["marker"]

        self.data[0]["marker"]["size"] = [marker_styles["default_marker_size"] for k in markers["size"]]
        self.data[0]["marker"]["color"] = [marker_styles["default_marker_color"] for k in markers["color"]]

        # if cd1 is not None:
        #     markers_size[cd1['points'][0]['pointNumber']] = self.clicked_marker_size
        #     markers_color[cd1['points'][0]['pointNumber']] = self.clicked_marker_color

        if hd is not None:
            # markers_size[hd['points'][0]['pointNumber']] = self.hovered_marker_size
            # markers_color[hd['points'][0]['pointNumber']] = self.hovered_marker_color
            self.data[0]["marker"]["size"][hd['points'][0]['pointNumber']] = marker_styles["hovered_marker_size"]
            self.data[0]["marker"]["color"][hd['points'][0]['pointNumber']] = marker_styles["hovered_marker_color"]

        # markers["size"] = tuple(markers_size)
        # markers["color"] = tuple(markers_color)
        # self.data[0]["marker"] = markers



class Markers:

    def __init__(self, marker_dict, x=None, y=None, z=None):
        self.size = marker_dict["size"]
        self.color = marker_dict["color"]
        self.x = x
        self.y = y
        self.z = z

    def update(self, hd, cd_own=None, cd_foreign=None, labels_own=None, labels_foreign=None, marker_styles=None):

        self.size = [marker_styles["default_marker_size"] for k in self.size]
        self.color = [marker_styles["default_marker_color"] for k in self.color]

        if cd_own is not None:
            self.size[cd_own['points'][0]['pointNumber']] = marker_styles["clicked_marker_size"]
            self.color[cd_own['points'][0]['pointNumber']] = marker_styles["clicked_marker_color"]

        if cd_foreign is not None:

            # get clickedData from the foreign chart
            x_foreign = cd_foreign["points"][0]["x"]
            y_foreign = cd_foreign["points"][0]["y"]
            z_foreign = cd_foreign["points"][0]["z"]

            # get axis names from labels_own and labels_foreign
            o_xl, o_yl, o_zl = labels_own[0], labels_own[1], labels_own[2]
            o_zerol = labels_own[3]
            f_xl, f_yl, f_zl = labels_foreign[0], labels_foreign[1], labels_foreign[2]

            # create dict with coordinates of the points from the own graph
            solutions_dict = {o_xl: self.x, o_yl: self.y, o_zl: self.z, o_zerol: [0 for i in self.x]}

            # create check list necessary to validate points
            if (f_xl == "mx") or (f_xl == "nx"):
                x_check = ["px", "sx"]
            elif (f_xl == "my") or (f_xl == "ny"):
                x_check = ["py", "sy"]
            elif (f_xl == "px") or (f_xl == "sx"):
                x_check = ["mx", "nx"]
            elif (f_xl == "py") or (f_xl == "sy"):
                x_check = ["my", "ny"]
            if (f_yl == "mx") or (f_yl == "nx"):
                y_check = ["px", "sx"]
            elif (f_yl == "my") or (f_yl == "ny"):
                y_check = ["py", "sy"]
            elif (f_yl == "px") or (f_yl == "sx"):
                y_check = ["mx", "nx"]
            elif (f_yl == "py") or (f_yl == "sy"):
                y_check = ["my", "ny"]
            if (f_zl == "mx") or (f_zl == "nx"):
                z_check = ["px", "sx"]
            elif (f_zl == "my") or (f_zl == "ny"):
                z_check = ["py", "sy"]
            elif (f_zl == "px") or (f_zl == "sx"):
                z_check = ["mx", "nx"]
            elif (f_zl == "py") or (f_zl == "sy"):
                z_check = ["my", "ny"]

            # check conditionals and append list of the rejected points
            rejected_list = []
            for i in range(len(self.x)):
                conditionals = [
                    f"{solutions_dict[x_check[0]][i]:.1f}" == f"{x_foreign:.1f}",
                    f"{solutions_dict[x_check[1]][i]:.1f}" == f"{x_foreign:.1f}",
                    f"{solutions_dict[y_check[0]][i]:.1f}" == f"{y_foreign:.1f}",
                    f"{solutions_dict[y_check[1]][i]:.1f}" == f"{y_foreign:.1f}",
                    f"{solutions_dict[z_check[0]][i]:.1f}" == f"{z_foreign:.1f}",
                    f"{solutions_dict[z_check[1]][i]:.1f}" == f"{z_foreign:.1f}"
                ]
                if any(conditionals):
                    rejected_list.append(i)

            # change styles the markers from the rejected list
            for index in rejected_list:
                self.color[index] = marker_styles["rejected_marker_color"]
                self.size[index] = marker_styles["rejected_marker_size"]

        if hd is not None:
            self.size[hd['points'][0]['pointNumber']] = marker_styles["hovered_marker_size"]
            self.color[hd['points'][0]['pointNumber']] = marker_styles["hovered_marker_color"]

        return dict(size=self.size, color=self.color, opacity=1,
                    line=dict(width=marker_styles["default_marker_edge_size"],
                              color=marker_styles["default_marker_edge_color"]))



