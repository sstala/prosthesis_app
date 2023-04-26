import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_table
from dash.dependencies import Input, Output, State
from dash import ctx, no_update
import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import decimal
# from plotly.callbacks import Points, InputDeviceState
#
# points, state = Points(), InputDeviceState()
# import plotly.io as pio
# svg_renderer = pio.renderers["svg"]
# svg_renderer.engine = 'kaleido'  # static image generation dependency, install it using pip or conda.
# pio.renderers.default = "svg"
# from abc import ABC, abstractmethod
from plotclass import Markers

marker_styles = dict(
    FindParamsScatter3d=dict(
        default_marker_size=4,
        default_marker_color="#85c5ed",
        clicked_marker_size=12,
        clicked_marker_color="#67cf1d",
        hovered_marker_size=10,
        hovered_marker_color="#1dcfb1",
        default_marker_edge_size=8,
        default_marker_edge_color='#1573ad',
        rejected_marker_size=6,
        rejected_marker_color="red",
    ),
    MechanismDimensions=dict(
        default_marker_size=8,
        default_marker_color="#a32c2c",
        clicked_marker_size=4,
        clicked_marker_color="black",
        hovered_marker_size=10,
        hovered_marker_color="#ed4242",
        default_marker_edge_size=1,
        default_marker_edge_color='#ffd64f',
        rejected_marker_size=4,
        rejected_marker_color="black",
    ),
)


def find_decimals(value):
    return (abs(decimal.Decimal(str(value)).as_tuple().exponent))


def init_global_vals():
    global table_header_names, licznik, \
        first_time, tabstyle_header_enabled, tabstyle_cell_enabled, tabstyle_header_disabled, tabstyle_cell_disabled
    licznik = 0
    table_header_names = ['mx', 'my', 'nx', 'ny', 'px', 'py', 'sx', 'sy']
    first_time = True
    tabstyle_header_enabled = {'color': "#809ead"}
    tabstyle_cell_enabled = {
        'text-align': 'center',
        'background-color': '#0a2836',
        'border': '0px',
        'color': '#ffd64f',
        'minWidth': '30px', 'width': '30px', 'maxWidth': '30px',
    }
    tabstyle_header_disabled = {'color': "#8e9091"}
    tabstyle_cell_disabled = {
        'text-align': 'center',
        'background-color': '#0a2836',
        'border': '0px',
        'color': '#8e9091',
        'minWidth': '30px', 'width': '30px', 'maxWidth': '30px',
        'cursor': 'not-allowed'
    }


def prepare_content_schema_in_tab1():
    global tabstyle_header_disabled, tabstyle_header_disabled, tabstyle_header_enabled, tabstyle_header_enabled
    return ([
        html.Div(
            id='other_tab1',
            children=[
                html.Div(
                    id='set-val-panel',
                    children=[
                        html.Div([
                            html.Label('Enter the segment lengths'),
                        ], className="set-val-panel-row-label"),
                        html.Div([
                            html.Label('|AB|'),
                            html.Label('='),
                            dcc.Input(
                                id="id_AB",
                                value=40,
                                placeholder='value',
                                type='number',
                                readOnly=False,
                                disabled=True,
                            ),
                            html.Label('25 < |AB| < 50 [mm]', style={'font-size': 10}),
                        ], className="set-val-panel-row"),
                        html.Div([
                            html.Label('|BC|'),
                            html.Label('='),
                            dcc.Input(
                                id="id_BC",
                                value=35,
                                placeholder='value',
                                type='number',
                                readOnly=False,
                                disabled=True,
                            ),
                            html.Label('20 < |BC| < 45 [mm]', style={'font-size': 10}),
                        ], className="set-val-panel-row"),
                        html.Div([
                            html.Label('|CD|'),
                            html.Label('='),
                            dcc.Input(
                                id="id_CD",
                                value=25,
                                placeholder='value',
                                type='number',
                                readOnly=False,
                                disabled=True,
                            ),
                            html.Label('10 < |CD| < 30 [mm]', style={'font-size': 10}),
                        ], className="set-val-panel-row"),
                        html.Div([
                            html.Label('Choose variables equal zero'),
                        ], className="set-val-panel-row-label"),
                        html.Div([
                            html.Div([
                                dcc.Dropdown(
                                    [
                                        {
                                            "label": "mx", "value": "mx",
                                            'disabled': False
                                        },
                                        {
                                            "label": "my", "value": "my",
                                            'disabled': True
                                        },
                                        {
                                            "label": "nx", "value": "nx",
                                            'disabled': True
                                        },
                                        {
                                            "label": "ny", "value": "ny",
                                            'disabled': True
                                        }
                                    ],
                                    id='id_dd1',
                                    searchable=False,
                                    disabled=True,
                                    value="mx",
                                    placeholder="Choose variable",
                                ),

                            ], className="dropdown-div-disabled", id="id_dropdown1_div", ),

                            html.Div([
                                dcc.Dropdown(
                                    [
                                        {
                                            "label": "px", "value": "px",
                                            'disabled': False
                                        },
                                        {
                                            "label": "py", "value": "py",
                                            'disabled': True
                                        },
                                        {
                                            "label": "sx", "value": "sx",
                                            'disabled': True
                                        },
                                        {
                                            "label": "sy", "value": "sy",
                                            'disabled': True
                                        }
                                    ],
                                    id='id_dd2',
                                    searchable=False,
                                    disabled=True,
                                    placeholder="Choose variable",
                                    value="px",
                                )
                            ], className="dropdown-div-disabled", id="id_dropdown2_div")
                        ], className="set-val-panel-row"),
                        html.Div([
                            dbc.Button(
                                'Calculate',
                                id='id_btn',
                                n_clicks=0,
                                disabled=True,
                                className='my-button',
                            )
                        ], className="set-val-panel-btn"),
                        html.Div([
                            html.Label('Select solutions from the graphs'),
                        ], className="set-val-panel-row-label"),
                        html.Div([
                            html.Div(
                                id='table_div',
                                children=[
                                    dash_table.DataTable(
                                        id='id_table_show_sols0',
                                        columns=([{'id': p, 'name': p} for p in table_header_names]),
                                        data=[dict(**{param: '-' for param in table_header_names})],
                                        style_header=tabstyle_header_disabled,
                                        style_cell=tabstyle_cell_disabled,
                                        style_data_conditional=[{
                                            "if": {"state": "selected"},
                                            "backgroundColor": "inherit !important",
                                            "border": "inherit !important",
                                        }]
                                    )
                                ],
                                className="table_show_sols"),

                        ], className="table_show_sols3"),
                        html.Div([
                            dbc.Button(
                                "SHOW",
                                id='show_btn',
                                n_clicks=0,
                                disabled=True,
                                className="my-button"

                            ),
                            dbc.Button(
                                "\u21A9 ",
                                id='return_btn',
                                disabled=True,
                                n_clicks=0,
                                className="my-button"
                            )

                        ], className="set-val-panel-btn2"),

                    ]
                ),
                html.Div(
                    id='cont1',
                    children=[
                        html.Div(
                            id='cont2',
                            children=[
                                dcc.Tabs(
                                    value="dimensional_drawing",
                                    id="cont-tabs",
                                    className="cont-tab-cont",
                                    children=[
                                        dcc.Tab(
                                            id="cont-tab1",
                                            label="dimensional drawing",
                                            value="dimensional_drawing",
                                            className="cont-tab",
                                            selected_className="cont-tab--selected",
                                            disabled=False,
                                        ),
                                        dcc.Tab(
                                            id="cont-tab2",
                                            label="choose solutions",
                                            value="choose_solutions",
                                            className="cont-tab",
                                            selected_className="cont-tab--selected",
                                            disabled=False,
                                        ),
                                        dcc.Tab(
                                            id="cont-tab3",
                                            label="model parameters",
                                            value="model_parameters",
                                            className="cont-tab",
                                            selected_className="cont-tab--selected",
                                            disabled=True,
                                            disabled_className="cont-tab--disabled"
                                        ),
                                    ]
                                )
                            ]
                        ),
                        html.Div(
                            id='cont3',
                            children=[
                                html.Div(
                                    id="id_pct_dim1",
                                    children=[
                                        html.Img(id="id_pct_dim2", src=app.get_asset_url("wymiary.png")),
                                    ]
                                ),
                            ]
                        )
                    ]
                )

            ]
        ),

    ]
    )


def build_main_tabs():
    return (
        html.Div(
            id="menu-tabs",
            children=[
                dcc.Tabs(
                    id="id_app-tabs",
                    parent_className='custom-tabs',
                    className="custom-tabs-container",
                    children=[
                        dcc.Tab(
                            id="id_tab0",
                            label="INFO",
                            value="tab0",
                            className="custom-tab",
                            selected_className="custom-tab--selected"
                        ),
                        dcc.Tab(
                            id="id_tab1",
                            label="TAB1",
                            value="tab1",
                            className="custom-tab",
                            selected_className="custom-tab--selected"
                        ),
                        dcc.Tab(
                            id="id_tab2",
                            label="TAB2",
                            value="tab2",
                            className="custom-tab",
                            selected_className="custom-tab--selected"
                        ),
                        dcc.Tab(
                            id="id_tab3",
                            label="TAB3",
                            value="tab3",
                            className="custom-tab",
                            selected_className="custom-tab--selected"
                        ),
                    ]
                )
            ]

        )
    )


def calculate_model(AB, BC, CD, table_dict):
    printinfo("calculate_model", "-")

    mx = table_dict["mx"]
    my = table_dict["my"]
    nx = table_dict["nx"]
    ny = table_dict["ny"]
    px = table_dict["px"]
    py = table_dict["py"]
    sx = table_dict["sx"]
    sy = table_dict["sy"]

    def auxliary_edge_length(p1x, p2x, p1y, p2y):
        return ((p1x - p2x) ** 2 + (p1y - p2y) ** 2) ** 0.5

    def triangle_perimeter(len1, len2, len3):
        return (len1 + len2 + len3) / 2

    def triangle_area(perimeter, len1, len2, len3):
        return (perimeter * (perimeter - len1) * (perimeter - len2) * (perimeter - len3)) ** 0.5

    def round_arcus(x):
        if x < -1:
            return -1
        elif x > 1:
            return 1
        else:
            return x

    # mx=0
    # my = 5.5
    # nx=5.5
    # ny=0
    # px=0
    # py=5.5
    # sx=3.1
    # sy=3.5
    # AB = 40
    # BC = 35
    # CD = 25

    Ax = 0
    Ay = 0
    AB = AB
    BC = BC
    CD = CD
    Ex = 0 + mx
    Ey = 0 - my
    Fx0 = AB - nx
    Fy0 = ny
    Bx0 = AB
    By0 = 0
    Gx0 = AB - px
    Gy0 = -py
    Hx0 = AB + BC - sx
    Hy0 = sy
    Cx0 = AB + BC
    Cy0 = 0
    Dx0 = AB + BC + CD
    Dy0 = 0
    BG = (px ** 2 + py ** 2) ** 0.5
    BF = (nx ** 2 + ny ** 2) ** 0.5
    EF = ((AB - mx - nx) ** 2 + (my + ny) ** 2) ** 0.5

    GH = ((py + sy) ** 2 + (BC - sx) ** 2) ** 0.5
    CH = (sx ** 2 + sy ** 2) ** 0.5
    AE = auxliary_edge_length(Ax, Ex, Ay, Ey)

    joint_trajectory_dict = {}
    joint_trajectory_dict_keys = ['kABC', 'kBCD', 'alpha', 'Ax', 'Ay', 'Bx', 'By', 'Cx', 'Cy', 'Dx', 'Dy', 'Ex', 'Ey',
                                  'Fx', 'Fy', 'Gx', 'Gy', 'Hx', 'Hy', 'dB', 'dBx', 'dBy', 'dC', 'dCx', 'dCy', 'dD',
                                  'dDx', 'dDy', 'dE', 'dF', 'dFx', 'dFy', 'dG', 'dGx', 'dGy', 'dH', 'dHx', 'dHy']
    for key in joint_trajectory_dict_keys:
        joint_trajectory_dict[key] = []

    new_row = {}
    # --- dodałem nowy kod do obliczeni CF i DH
    CF = auxliary_edge_length(Cx0, Fx0, Cy0, Fy0)  #
    perCF = triangle_perimeter(BF, BC, CF)  #
    areCF = triangle_area(perCF, BF, BC, CF)  #
    gamCF = 180 - np.arcsin(areCF / (0.5 * BF * BC)) * (180 / np.pi)  #
    DH = auxliary_edge_length(Dx0, Hx0, Dy0, Hy0)
    perDH = triangle_perimeter(CH, DH, CD)  #
    areDH = triangle_area(perDH, CH, DH, CD)  #
    gamDH = 180 - np.arcsin(areDH / (0.5 * CD * CH)) * (180 / np.pi)  #

    step = 0.1
    for i in np.arange(0, 90.001, step):
        alfa = i
        Bx = Ax + AB * np.cos(alfa * np.pi / 180)
        By = Ay - AB * np.sin(alfa * np.pi / 180)

        BE = auxliary_edge_length(Ex, Bx, Ey, By)

        per0 = triangle_perimeter(AE, AB, BE)
        are0 = triangle_area(per0, AE, AB, BE)
        gam0 = 90 + alfa - np.arcsin(are0 / (0.5 * AB * BE)) * (180 / np.pi)

        per1 = triangle_perimeter(BF, EF, BE)
        are1 = triangle_area(per1, BF, EF, BE)
        gam1 = np.arcsin(are1 / (0.5 * EF * BE)) * (180 / np.pi)

        Fx = Ex + EF * np.sin((gam0 - gam1) * np.pi / 180)
        Fy = Ey + EF * np.cos((gam0 - gam1) * np.pi / 180)

        AF = auxliary_edge_length(Ax, Fx, Ay, Fy)
        # perABF = triangle_perimeter(AB, AF, BF)
        # areABF = triangle_area(perABF, AB, AF, BF)
        # kAFB = np.arcsin(areABF / (0.5 * AB * BF)) * (180 / np.pi)
        # kAFB2 = 90 - np.arccos(areABF / (0.5 * AB * BF)) * (180 / np.pi)

        # costam = round((AF**2 - BF**2 - AB**2)/(-2*BF*AB),4)
        round_arcus

        # kAFB3 = np.arccos(round((AF ** 2 - BF ** 2 - AB ** 2) / (-2 * BF * AB), 14)) * (180 / np.pi)
        kAFB3 = np.arccos(round_arcus((AF ** 2 - BF ** 2 - AB ** 2) / (-2 * BF * AB))) * (180 / np.pi)  # tu

        tet1 = kAFB3 + alfa  # 2

        Gx = Bx - BG * np.sin(alfa * np.pi / 180)
        Gy = By - BG * np.cos(alfa * np.pi / 180)

        # stara definicja punktu C
        # Cx = Bx + BC * np.cos(tet1 * np.pi / 180)
        # Cy = By - BC * np.cos((90 - tet1) * np.pi / 180)

        # nowa definicja punktu C
        Cx = Bx + BC * np.cos((tet1 + gamCF - 180) * np.pi / 180)
        Cy = By + BC * np.cos((90 + tet1 + gamCF - 180) * np.pi / 180)
        # print(f'a: {i}   x:{tet1+gamCF-180}   cosx:{np.cos((tet1+gamCF-180) * np.pi/180)}   y:{tet1+gamCF-90}   cosy:{np.cos((tet1+gamCF-90) * np.pi/180)}')
        # print(f'AF: {AF}   BF: {BF}   costam: {costam}   kAFB3: {kAFB3}   alfa: {i}   tet1: {tet1}')

        # tet2 = np.arccos(np.round_(((Cx - Bx) / BC), 4)) * 180 / np.pi

        CG = auxliary_edge_length(Cx, Gx, Cy, Gy)
        per2 = triangle_perimeter(BG, BC, CG)
        are2 = triangle_area(per2, BG, BC, CG)
        # gam21 = np.arcsin(round(are2 / (0.5 * BC * CG), 14)) * 180 / np.pi
        # gam22 = np.arcsin(round(are2 / (0.5 * BG * BC), 14)) * 180 / np.pi
        gam21 = np.arcsin(round_arcus(are2 / (0.5 * BC * CG))) * 180 / np.pi  # tu
        gam22 = np.arcsin(round_arcus(are2 / (0.5 * BG * BC))) * 180 / np.pi  # tu
        gam2 = 180 - gam21 - gam22

        # GH = ((5.5 + 3.5) ** 2 + (35 - 3.1) ** 2) ** 0.5
        # CH = (3.1 ** 2 + 3.5 ** 2) ** 0.5
        # GH = ((py + sy)**2 + (BC - sx)**2)**0.5
        # CH = (sx ** 2 + sy ** 2) ** 0.5

        per3 = triangle_perimeter(GH, CH, CG)
        are3 = triangle_area(per3, GH, CH, CG)
        gam3 = np.arcsin(are3 / (0.5 * GH * CG)) * 180 / np.pi

        tet3 = gam2 - gam3 + alfa

        Hx = Gx + GH * np.sin(tet3 * np.pi / 180)
        Hy = Gy + GH * np.cos(tet3 * np.pi / 180)

        BH = auxliary_edge_length(Bx, Hx, By, Hy)

        # per4 = triangle_perimeter(BH, CH, BC)
        # are4 = triangle_area(per4, BH, CH, BC)
        gam4 = np.arccos((CH ** 2 + BC ** 2 - BH ** 2) / (2 * CH * BC)) * 180 / np.pi

        DH = ((CD + sx) ** 2 + sy ** 2) ** 0.5

        # per5 = triangle_perimeter(CH, CD, DH)
        # are5 = triangle_area(per5, CH, CD, DH)
        # gam5 = np.arccos(round((CH ** 2 + CD ** 2 - DH ** 2) / (2 * CH * CD), 14)) * 180 / np.pi
        gam5 = np.arccos(round_arcus((CH ** 2 + CD ** 2 - DH ** 2) / (2 * CH * CD))) * 180 / np.pi

        tet4 = gam5 + gam4 - (90 - (tet1 - 180 + gamCF))

        # stara definicja punktu D
        # Dx = Cx + CD * np.sin(tet4 * np.pi / 180)
        # Dy = Cy + CD * np.cos(tet4 * np.pi / 180)

        # stara definicja punktu D
        # tet5 = np.arccos(np.round_(((Cx - Hx) / CH), 4)) * 180 / np.pi

        tet6 = -180 + gamCF + tet1 + gamDH + tet4
        tet7 = 360 - gam4 - gamDH - tet1 + 180 - gamCF
        tet8 = 180 - gamDH - gam4 + tet1
        # Dx = Cx + CD * np.cos((tet5 + gamDH - 180) * np.pi / 180)
        # Dy = Cy + CD * np.cos((90 + tet5 + gamDH - 180) * np.pi / 180)
        Dx = Cx + CD * np.sin(tet4 * np.pi / 180)
        Dy = Cy + CD * np.cos(tet4 * np.pi / 180)
        # print(f'gam21: {gam21}   gam22: {gam22}   BH: {BH} CH: {CH}   BC: {BC}   my: {my}    gam5: {gam5}   gam4: {gam4}    tet1: {(CH ** 2 + BC ** 2 - BH ** 2)}   CD:{(2 * CH * BC)}   tet4:{tet4}   sin:{np.sin(tet4 * np.pi / 180)}   cos:{np.cos(tet4 * np.pi / 180)} ')

        AC = auxliary_edge_length(Ax, Cx, Ay, Cy)
        kABC = np.arccos(round_arcus((AC ** 2 - AB ** 2 - BC ** 2) / (-2 * AB * BC))) * 180 / np.pi
        BD = auxliary_edge_length(Bx, Dx, By, Dy)
        kBCD = np.arccos(round_arcus((BD ** 2 - BC ** 2 - CD ** 2) / (-2 * BC * CD))) * 180 / np.pi
        if new_row:
            old_row = new_row
        else:
            old_row = {'kABC': kABC, 'kBCD': kBCD, 'alpha': i, 'Ax': Ax, 'Ay': Ay, 'Bx': Bx, 'By': By,
                       'Cx': Cx, 'Cy': Cy, 'Dx': Dx, 'Dy': Dy, 'Ex': Ex, 'Ey': Ey, 'Fx': Fx, 'Fy': Fy,
                       'Gx': Gx, 'Gy': Gy, 'Hx': Hx, 'Hy': Hy, }

        new_row = {'kABC': kABC, 'kBCD': kBCD, 'alpha': i, 'Ax': Ax, 'Ay': Ay, 'Bx': Bx, 'By': By, 'Cx': Cx, 'Cy': Cy,
                   'Dx': Dx, 'Dy': Dy, 'Ex': Ex, 'Ey': Ey, 'Fx': Fx, 'Fy': Fy, 'Gx': Gx, 'Gy': Gy, 'Hx': Hx, 'Hy': Hy,
                   'dB': round(((((old_row['Bx'] - Bx) ** 2 + (old_row['By'] - By) ** 2) ** 0.5) / step), 15),
                   'dC': round(((((old_row['Cx'] - Cx) ** 2 + (old_row['Cy'] - Cy) ** 2) ** 0.5) / step), 15),
                   'dD': round(((((old_row['Dx'] - Dx) ** 2 + (old_row['Dy'] - Dy) ** 2) ** 0.5) / step), 15),
                   'dE': round(((((old_row['Ex'] - Ex) ** 2 + (old_row['Ey'] - Ey) ** 2) ** 0.5) / step), 15),
                   'dF': round(((((old_row['Fx'] - Fx) ** 2 + (old_row['Fy'] - Fy) ** 2) ** 0.5) / step), 15),
                   'dG': round(((((old_row['Gx'] - Gx) ** 2 + (old_row['Gy'] - Gy) ** 2) ** 0.5) / step), 15),
                   'dH': round(((((old_row['Hx'] - Hx) ** 2 + (old_row['Hy'] - Hy) ** 2) ** 0.5) / step), 15),
                   'dBx': round(((Bx - old_row['Bx']) / step), 15), 'dBy': round(((By - old_row['By']) / step), 15),
                   'dCx': round(((Cx - old_row['Cx']) / step), 15), 'dCy': round(((Cy - old_row['Cy']) / step), 15),
                   'dDx': round(((Dx - old_row['Dx']) / step), 15), 'dDy': round(((Dy - old_row['Dy']) / step), 15),
                   'dFx': round(((Fx - old_row['Fx']) / step), 15), 'dFy': round(((Fy - old_row['Fy']) / step), 15),
                   'dGx': round(((Gx - old_row['Gx']) / step), 15), 'dGy': round(((Gy - old_row['Gy']) / step), 15),
                   'dHx': round(((Hx - old_row['Hx']) / step), 15), 'dHy': round(((Hy - old_row['Hy']) / step), 15),
                   }
        for key in joint_trajectory_dict_keys:
            joint_trajectory_dict[key].append(new_row[key])

    max_abs_value = 0
    for key in ["dBx", "dBy", "dCx", "dCy", "dDx", "dDy", "dFx", "dFy", "dGx", "dGy", "dHx", "dHy"]:
        if max_abs_value < np.max(np.abs(joint_trajectory_dict[key])):
            max_abs_value = np.max(np.abs(joint_trajectory_dict[key]))

    for key in ["Bx", "Cx", "Dx", "Fx", "Gx", "Hx"]:
        joint_trajectory_dict["vh" + key] = [
            joint_trajectory_dict[key][i] + 40 * joint_trajectory_dict["d" + key][i] / max_abs_value
            for i in range(len(joint_trajectory_dict[key]))]
        joint_trajectory_dict["vv" + key] = [
            joint_trajectory_dict[key][i] for i in range(len(joint_trajectory_dict[key]))]

        joint_trajectory_dict["vhs" + key[0]] = [
            ["circle", "triangle-left"] if element < 0
            else ["circle", "triangle-right"] if element > 0
            else ["circle", "circle"]
            for element in joint_trajectory_dict["d" + key]
        ]

    for key in ["By", "Cy", "Dy", "Fy", "Gy", "Hy"]:
        joint_trajectory_dict["vh" + key] = [
            joint_trajectory_dict[key][i] for i in range(len(joint_trajectory_dict[key]))]
        joint_trajectory_dict["vv" + key] = [
            joint_trajectory_dict[key][i] + 40 * joint_trajectory_dict["d" + key][i] / max_abs_value
            for i in range(len(joint_trajectory_dict[key]))]

        joint_trajectory_dict["vvs" + key[0]] = [
            ["circle", "triangle-down"] if element < 0
            else ["circle", "triangle-up"] if element > 0
            else ["circle", "circle"]
            for element in joint_trajectory_dict["d" + key]
        ]

    for key in ["dB", "dBx", "dBy", "dC", "dCx", "dCy", "dD", "dDx", "dDy", "dF", "dFx", "dFy", "dG", "dGx", "dGy",
                "dH", "dHx", "dHy"]:
        joint_trajectory_dict["tv" + key] = [f"{np.abs(element):.3f}" for element in joint_trajectory_dict[key]]
        joint_trajectory_dict["mv" + key] = [np.abs(element) / max_abs_value for element in joint_trajectory_dict[key]]

    return joint_trajectory_dict


def build_find_parameters_plot(ab_length, bc_length, cd_length, eq_zero_var1, eq_zero_var2):
    printinfo("build_find_parameters_plot", "-")

    def prepare_data():

        def calculate_solutions():
            math_model = dict(
                mx="(-ab * (e + i + j) + 2 * e * i) / (-ab - 2 * j)",
                my="(-ab * (i + j - e) + 2 * e * j) / (ab - 2 * i)",
                nx="(-ab * (i + e - j) + 2 * i * e) / (ab - 2 * j)",
                ny="(-ab * (i + e - j) + 2 * j * ny) / (ab - 2 * i)",
                px="(bc * (e - i - j) - 2 * i * e) / (-bc - 2 * j)",
                py="(bc * (-i - j + e) + 2 * j * e) / (-bc + 2 * i)",
                sx="(bc * (i - e + j) - 2 * e * i) / (bc - 2 * j)",
                sy="(bc * (i - e + j) + 2 * e * j) / (bc + 2 * i)",
            )
            labels = dict(
                mx=["nx", "ny", "my"],
                my=["nx", "ny", "mx"],
                nx=["mx", "my", "ny"],
                ny=["mx", "my", "nx"],
                px=["sx", "sy", "py"],
                py=["sx", "sy", "px"],
                sx=["px", "py", "sy"],
                sy=["px", "py", "sx"],
            )

            def calculate(eq_zero_var):
                stp, min_val, max_val = 0.1, 0, 8

                ab, bc, cd = float(ab_length), float(bc_length), float(cd_length)

                e = 0
                x = np.arange(min_val, max_val, stp)
                y = np.arange(min_val, max_val, stp)
                x_filtered = []
                y_filtered = []
                z_filtered = []
                for i in x:
                    for j in y:
                        z = eval(math_model[eq_zero_var])
                        if find_decimals(z) == 1:
                            x_filtered.append(i)
                            y_filtered.append(j)
                            z_filtered.append(z)
                return x_filtered, labels[eq_zero_var][0], y_filtered, labels[eq_zero_var][1], \
                    z_filtered, labels[eq_zero_var][2]

            x1, x1_label, y1, y1_label, z1, z1_label = calculate(eq_zero_var1)
            x2, x2_label, y2, y2_label, z2, z2_label = calculate(eq_zero_var2)

            return x1, x1_label, y1, y1_label, z1, z1_label, x2, x2_label, y2, y2_label, z2, z2_label

        x1, x1_l, y1, y1_l, z1, z1_l, x2, x2_l, y2, y2_l, z2, z2_l = calculate_solutions()

        plot1 = dict(x=x1, y=y1, z=z1, x_label=x1_l, y_label=y1_l, z_label=z1_l)
        plot2 = dict(x=x2, y=y2, z=z2, x_label=x2_l, y_label=y2_l, z_label=z2_l)

        return plot1, plot2

    def prepare_dict(x, y, z, x_label, y_label, z_label):
        """ FUNCTION FLOW:
        1-1. prepare_data_and_layout() -> prepare_hover_tags() -> prepare_axis_names()// get axis names
        1-2. prepare_data_and_layout() -> prepare_hover_tags() // get hovertext
        2-1. prepare_data_and_layout() -> prepare_scene_properties() -> prepare_axis_names()// get axis names
        2-2. prepare_data_and_layout() -> prepare_scene_properties() // get scene properties
        """

        def prepare_axis_names():
            xl = x_label[0] + '<sub>' + x_label[1] + '</sub>'
            yl = y_label[0] + '<sub>' + y_label[1] + '</sub>'
            zl = z_label[0] + '<sub>' + z_label[1] + '</sub>'
            return xl, yl, zl

        def prepare_hover_tags():
            """Function to prepare hovertext"""
            tags = []
            xl, yl, zl = prepare_axis_names()
            for i in range(len(x)):
                tag = '<b><i>' + xl + ':</b></i> ' + ' ' * (4 - len("{0:.1f}".format(x[i]))) + \
                      "{0:.1f}".format(x[i]) + ' [mm]<br>' + \
                      '<b><i>' + yl + ':</b></i> ' + ' ' * (4 - len("{0:.1f}".format(y[i]))) + \
                      "{0:.1f}".format(y[i]) + ' [mm]<br>' + \
                      '<b><i>' + zl + ':</b></i> ' + ' ' * (4 - len("{0:.1f}".format(z[i]))) + \
                      "{0:.1f}".format(z[i]) + ' [mm]' + '<extra></extra>'
                tags.append(tag)
            return tags

        def prepare_scene_properties():
            scene_dict = {}
            axis_names = list(prepare_axis_names())
            for idx, name in enumerate(["xaxis", "yaxis", "zaxis"]):
                scene_dict[name] = dict(nticks=8, title=axis_names[idx], backgroundcolor='rgba(0, 0, 0, 0)',
                                        gridcolor='white', gridwidth=1, color='white')
            scene_dict["camera"] = dict(up={'x': 0, 'y': 0, 'z': 1}, center={'x': 0, 'y': 0, 'z': -0.1},
                                        eye={'x': -1, 'y': -1, 'z': 1})
            scene_dict["aspectratio"] = {"x": 0.7, "y": 0.7, "z": 0.6}
            return scene_dict

        def prepare_data_and_layout():
            data = [go.Scatter3d(
                x=x, y=y, z=z, type='scatter3d', mode="markers", hovertemplate=prepare_hover_tags(),
                hoverlabel=dict(bgcolor='#103d52', bordercolor='#ffd64f', font=dict(family="Consolas")),
                marker=dict(opacity=1,
                            size=[marker_styles["FindParamsScatter3d"]["default_marker_size"]
                                  for k in range(len(x))],
                            color=[marker_styles["FindParamsScatter3d"]["default_marker_color"]
                                   for k in range(len(x))],
                            line=dict(width=marker_styles["FindParamsScatter3d"]["default_marker_edge_size"],
                                      color=marker_styles["FindParamsScatter3d"]["default_marker_edge_color"])))]
            layout = dict(margin={'l': 0, 'r': 0, 't': 0, 'b': 0}, plot_bgcolor='#111111',
                          paper_bgcolor='rgba(0,0,0,0)', scene=prepare_scene_properties())
            return data, layout

        return prepare_data_and_layout()

    plot1, plot2 = prepare_data()
    data1, layout1 = prepare_dict(plot1["x"], plot1["y"], plot1["z"], plot1["x_label"], plot1["y_label"],
                                  plot1["z_label"])
    data2, layout2 = prepare_dict(plot2["x"], plot2["y"], plot2["z"], plot2["x_label"], plot2["y_label"],
                                  plot2["z_label"])
    fig1 = go.Figure(data1, layout1)
    fig2 = go.Figure(data2, layout2)

    # return fig1, fig2
    return (dcc.Graph(clear_on_unhover=True, config={'displayModeBar': False, 'scrollZoom': False}, id='id_graph1',
                      figure=fig1),
            dcc.Graph(clear_on_unhover=True, config={'displayModeBar': False, 'scrollZoom': False}, id='id_graph2',
                      figure=fig2))


def build_mechanism_dimensions_chart(AB, BC, CD, eq_zero_var1, eq_zero_var2, table):
    from operator import itemgetter

    def prepare_hover_tags(point_list, x_list, y_list, l_tag):
        tags = []
        # if names in 'point_list' exist in 'long_tag' then name have long tag, otherwise short tag
        for idx, elm in enumerate(point_list):
            if elm not in l_tag:
                tagtxt = "<span style='font-size:25; color:#16b9de'>" + "<b>" + elm + "</b></span><br>" + \
                         "<b><i>x</i>:  </b>" + "<span style='color: #ed4242'>-</span>" + "  [mm]<br>" + \
                         "<b><i>Y</i>:  </b>" + "<span style='color: #ed4242'>-</span>" + "  [mm]<br>" + \
                         "<extra></extra>"
            else:
                xlen, ylen = len(x_list[idx]), len(y_list[idx])

                if xlen > ylen:
                    xl = "<b><i>x</i>:  </b>" + x_list[idx] + "  [mm]"
                    yl = "<b><i>y</i>:  </b>" + " " * (xlen - ylen) + y_list[idx] + "  [mm]"
                else:
                    xl = "<b><i>x</i>:  </b>" + " " * (ylen - xlen) + x_list[idx] + "  [mm]"
                    yl = "<b><i>y</i>:  </b>" + y_list[idx] + "  [mm]"

                tagtxt = "<span style='font-size:25; color:#16b9de'>" + "<b>" + elm + "</b></span><br>" + xl + \
                         "<br>" + yl + "<extra></extra>"
            tags.append(tagtxt)
        return tags

    def set_basic_dimensional_annotation():
        an_abc = ["<i><b>m<sub>y</sub></b></i>", "<i><b>n<sub>y</sub></b></i>", "<i><b>n<sub>x</sub></b></i>"]
        an_bcd = ["<i><b>p<sub>y</sub></b></i>", "<i><b>s<sub>y</sub></b></i>", "<i><b>s<sub>x</sub></b></i>"]
        return an_abc, an_bcd

    def set_dedicated_dimensional_annotation(tab, variables_names_list):
        annotation = []
        for elm in variables_names_list:
            if tab[0][elm] == "-":
                text = "<i><b>" + elm[0] + "<sub>" + elm[1] + "</sub></b></i>"
            else:
                text = "<i><b>" + elm[0] + "<sub>" + elm[1] + "</sub></b></i>" + \
                       " = " + "{:.1f}".format(tab[0][elm]) + "[mm]"
            annotation.append(text)
        return annotation

    def set_layout(ancl, x_range, annotation):
        layout = go.Layout(
            annotations=[
                go.layout.Annotation(x=ancl[0][0], y=ancl[0][1], xref="x", yref="y", text=annotation[0],
                                     align='left', valign="bottom",
                                     font=dict(color="#ffd64f", family='consolas', size=10),
                                     showarrow=False, yanchor='top', xanchor='right', textangle=-90),
                go.layout.Annotation(x=ancl[1][0], y=ancl[1][1], xref="x", yref="y", text=annotation[1],
                                     align='left', valign="bottom",
                                     font=dict(color="#ffd64f", family='consolas', size=10),
                                     showarrow=False, yanchor='top', xanchor='right', textangle=-90),
                go.layout.Annotation(x=ancl[2][0], y=ancl[2][1], xref="x", yref="y", text=annotation[2],
                                     align='left', valign="top",
                                     font=dict(color="#ffd64f", family='consolas', size=10),
                                     showarrow=False, yanchor='bottom', xanchor='right')
            ],
            hovermode=None, hoverlabel=None, showlegend=False, margin=dict(l=10, r=10, t=10, b=10),
            paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)',
            xaxis=dict(zeroline=False, gridwidth=1, gridcolor='#103d52', color='white', range=x_range),
            yaxis=dict(zeroline=False, scaleratio=1, scaleanchor="x", gridcolor='#103d52', gridwidth=1,
                       range=[-28, 28], color='white'),
        )
        return layout

    def prepare_dimensional_annotation(tab, l_tag, eq_zero_var1, eq_zero_var2, x_list1, x_list2, y_list1, y_list2):

        # Prepare joints names list
        point_list1 = ["A", "B", "C", "E", "F", "G"]
        point_list2 = ["F", "B", "C", "H", "D", "G"]
        # Prepare dimensional variable names when chart is showed first time
        if tab is None:
            annotation_abc, annotation_bcd = set_basic_dimensional_annotation()
            abc_tags, bcd_tags = [], []
        # Prepare dimensional variable names when points are hovered or clicked (only if chart is not showed first time)
        else:
            # the for loop appends to the list the joint names that are related to the parameters that exist in the
            # table
            for key, value in tab[0].items():
                if value != "-":
                    if key in ["mx", "my"]:
                        l_tag.append("E")
                    elif key in ["nx", "ny"]:
                        l_tag.append("F")
                    elif key in ["px", "py"]:
                        l_tag.append("G")
                    elif key in ["sx", "sy"]:
                        l_tag.append("H")

            abc_tags = prepare_hover_tags(point_list1, x_list1, y_list1, l_tag)
            bcd_tags = prepare_hover_tags(point_list2, x_list2, y_list2, l_tag)

            # prepare list with variables names // dimensions names
            variables_names_list1 = ["mx", "my", "ny", "nx"]
            variables_names_list1.remove(eq_zero_var1)
            variables_names_list2 = ["px", "py", "sy", "sx"]
            variables_names_list2.remove(eq_zero_var2)

            annotation_abc = set_dedicated_dimensional_annotation(table, variables_names_list1)
            annotation_bcd = set_dedicated_dimensional_annotation(table, variables_names_list2)
        return annotation_abc, annotation_bcd, abc_tags, bcd_tags

    def set_data(dim_dash_lines, dim_solid_lines, dim_arrows, segments_data, segment_lines, joints_markers, tags):
        # dimensional dash lines
        data1 = []
        for idx, elm in enumerate(dim_dash_lines):
            trace = go.Scatter(x=elm[0], y=elm[1], mode='lines', line=dict(color='#ffd64f', dash='dash', width=1),
                               hovertemplate=None, hoverinfo='skip')
            data1.append(trace)

        # dimensional solid lines
        data2 = []
        for idx, elm in enumerate(dim_solid_lines):
            trace = go.Scatter(x=elm[0], y=elm[1], mode='lines+markers',
                               line=dict(color='#ffd64f', dash='solid', width=1),
                               marker=dict(opacity=1, color='#ffd64f', line=dict(color='#ffd64f'),
                                           symbol=dim_arrows[idx], size=elm[2]),
                               hovertemplate=None, hoverinfo='skip')
            data2.append(trace)

        # segments
        data3 = []
        for idx, elm in enumerate(segments_data):
            trace = go.Scatter(x=elm[0], y=elm[1], mode='lines',
                               line=dict(color=segment_lines[idx], dash='solid', width=elm[2]),
                               hovertemplate=None, hoverinfo='skip')
            data3.append(trace)

        # joints
        data4 = [go.Scatter(x=joints_markers[0], y=joints_markers[1], mode="markers",
                            marker=dict(
                                opacity=1,
                                size=[marker_styles["MechanismDimensions"]["default_marker_size"] for k in
                                      joints_markers[0]],
                                color=[marker_styles["MechanismDimensions"]["default_marker_color"] for k in
                                       joints_markers[0]],
                                line=dict(
                                    color=marker_styles["MechanismDimensions"]["default_marker_edge_color"],
                                    width=marker_styles["MechanismDimensions"]["default_marker_edge_size"])),
                            hoverlabel=dict(
                                bgcolor='#103d52', bordercolor='#ffd64f', align="right", font=dict(family="Consolas")),
                            hovertemplate=tags, hoveron="points")]
        return data1 + data2 + data3 + data4

    # Define default parameters to show mechanism dimensions chart when any of points is not hovered or clicked.
    # Two elements in the dict are replaced by 0 depending on variables names equal zero
    mechanism_dimensions_chart_default_values = {"mx": 4, "my": 5, "nx": 3, "ny": 4, "px": 4, "py": 5, "sx": 2, "sy": 3,
                                                 eq_zero_var1: 0, eq_zero_var2: 0}

    # If table not exist then init dimensions variables by default values. Otherwise, use values from the table.
    # If variable is not exist then use default value.
    if table is None:
        mx, my, nx, ny, px, py, sx, sy = itemgetter("mx", "my", "nx", "ny", "px", "py", "sx", "sy") \
            (mechanism_dimensions_chart_default_values)
    else:
        new_dict = {}
        for els in ["mx", "my", "nx", "ny", "px", "py", "sx", "sy"]:
            if table[0][els] == "-":
                new_dict[els] = mechanism_dimensions_chart_default_values[els]
            else:
                new_dict[els] = table[0][els]
        mx, my, nx, ny, px, py, sx, sy = itemgetter("mx", "my", "nx", "ny", "px", "py", "sx", "sy")(new_dict)

    # Calculate the coordinates of the mechanism elements
    # calculate joints coordinates
    Ax, Ay = 0, 0
    Bx, By = Ax + AB, Ay
    Cx, Cy = Bx + BC, Ay
    Dx, Dy = Cx + CD, Ay
    Ex, Ey = Ax + mx, Ay - my
    Fx, Fy = Bx - nx, By + ny
    Gx, Gy = Bx - px, By - py
    Hx, Hy = Cx - sx, Cy + sy

    # calculate dimensions lines coordinates
    myd1x, myd1y = Ax - 10, Ey
    myd2x, myd2y = Ax - 10, Ay
    nxd1x, nxd1y = Fx, By - 15
    nxd2x, nxd2y = Bx, By - 15
    nyd1x, nyd1y = Cx + 10, Fy
    nyd2x, nyd2y = Cx + 10, Cy
    pyd1x, pyd1y = Bx - 10, By
    pyd2x, pyd2y = Bx - 10, Gy
    sxd1x, sxd1y = Hx, Cy - 15
    sxd2x, sxd2y = Cx, Cy - 15
    syd1x, syd1y = Dx + 10, Hy
    syd2x, syd2y = Dx + 10, Cy

    # prepare joints coordinates lists
    x_list1 = ["{0:.1f}".format(elm) for elm in [Ax, Bx, Cx, Ex, Fx, Gx]]
    y_list1 = ["{0:.1f}".format(elm) for elm in [Ay, By, Cy, Ey, Fy, Gy]]
    x_list2 = ["{0:.1f}".format(elm) for elm in [Fx, Bx, Cx, Hx, Dx, Gx]]
    y_list2 = ["{0:.1f}".format(elm) for elm in [Fy, By, Cy, Hy, Dy, Gy]]

    # list of base points, which are always defined according to the parameters entered in the input panel
    long_tag = ["A", "B", "C", "D"]

    annotation_abc, annotation_bcd, abc_tags, bcd_tags = prepare_dimensional_annotation(
        table, long_tag, eq_zero_var1, eq_zero_var2, x_list1, x_list2, y_list1, y_list2)

    # prepare parameters to build chart
    dim_dash_lines_abc = [[[Ex, myd1x], [Ey, myd1y]], [[Ax, myd2x], [Ay, myd2y]], [[Fx, nxd1x], [Fy, nxd1y]],
                          [[Bx, nxd2x], [By, nxd2y]], [[Fx, nyd1x], [Fy, nyd1y]], [[Bx, nyd2x], [By, nyd2y]]]

    dim_solid_lines_abc = [[[myd1x, myd2x, myd2x], [myd1y, myd2y, myd2y + 20], [7, 7, 0]],
                           [[nxd1x, nxd2x, nxd2x + 20], [nxd1y, nxd2y, nxd2y], [7, 7, 0]],
                           [[nyd1x, nyd1x, nyd2x], [nyd1y + 20, nyd1y, nyd2y], [0, 7, 7]]]

    dim_arrows_abc = [["arrow-up", "arrow-down", "arrow-down"],
                      ["arrow-right", "arrow-left", "arrow-down"],
                      ["arrow", "arrow-down", "arrow-up"]]

    segments_data_abc = [[[Ax, Bx, Gx], [Ay, By, Gy], 4], [[Fx, Bx, Cx], [Fy, By, Cy], 4], [[Ex, Fx], [Ey, Fy], 2]]

    segments_line_abc = ['blue', 'green', 'purple']

    joints_markers_abc = [[Ax, Bx, Cx, Ex, Fx, Gx], [Ay, By, Cy, Ey, Fy, Gy]]

    annotation_coordinates_list_abc = [[myd2x, myd2y + 20], [nyd1x, nyd1y + 20], [nxd2x + 20, nxd2y]]
    x_axis_range_abc = [-25, (AB + BC + 25)]

    mechanism_dimensions_chart_abc = go.Figure(
        data=set_data(dim_dash_lines_abc, dim_solid_lines_abc, dim_arrows_abc, segments_data_abc, segments_line_abc,
                      joints_markers_abc, abc_tags),
        layout=set_layout(annotation_coordinates_list_abc, x_axis_range_abc, annotation_abc))

    dim_dash_lines_bcd = [[[Bx, pyd1x], [By, pyd1y]], [[Gx, pyd2x], [Gy, pyd2y]], [[Hx, sxd1x], [Hy, sxd1y]],
                          [[Cx, sxd2x], [Cy, sxd2y]], [[Hx, syd1x], [Hy, syd1y]], [[Cx, syd2x], [Cy, syd2y]]]

    dim_solid_lines_bcd = [[[pyd2x, pyd1x, pyd2x], [pyd2y, pyd1y, pyd1y + 20], [7, 7, 0]],
                           [[sxd1x, sxd2x, sxd2x + 20], [sxd1y, sxd2y, sxd2y], [7, 7, 0]],
                           [[syd1x, syd1x, syd2x], [syd1y + 20, syd1y, syd2y], [0, 7, 7]]]

    dim_arrows_bcd = [["arrow-up", "arrow-down", "arrow-down"],
                      ["arrow-right", "arrow-left", "arrow-down"],
                      ["arrow", "arrow-down", "arrow-up"]]

    segments_data_bcd = [[[Fx, Bx, Cx], [Fy, By, Cy], 4], [[Hx, Cx, Dx], [Hy, Cy, Dy], 4], [[Gx, Hx], [Gy, Hy], 2]]

    segments_line_bcd = ['green', 'orange', 'purple']

    joints_markers_bcd = [[Fx, Bx, Cx, Hx, Dx, Gx], [Fy, By, Cy, Hy, Dy, Gy]]

    annotation_coordinates_list_bcd = [[pyd2x, pyd1y + 20], [syd1x, syd1y + 20], [sxd2x + 20, sxd2y]]

    x_axis_range_bcd = [(AB - 25), (AB + BC + CD + 25)]

    mechanism_dimensions_chart_bcd = go.Figure(
        data=set_data(dim_dash_lines_bcd, dim_solid_lines_bcd, dim_arrows_bcd, segments_data_bcd, segments_line_bcd,
                      joints_markers_bcd, bcd_tags),
        layout=set_layout(annotation_coordinates_list_bcd, x_axis_range_bcd, annotation_bcd))

    return dcc.Graph(clear_on_unhover=True,
                     config={'displayModeBar': False, 'scrollZoom': False},
                     id="id_graph7",
                     figure=mechanism_dimensions_chart_abc),\
        dcc.Graph(clear_on_unhover=True,
                  config={'displayModeBar': False, 'scrollZoom': False},
                  id="id_graph8",
                  figure=mechanism_dimensions_chart_bcd)

    # if table is None:
    #     mechanism_dimensions_chart_abc.data = []
    #     mechanism_dimensions_chart_abc.layout.annotations = []
    #     mechanism_dimensions_chart_bcd.data = []
    #     mechanism_dimensions_chart_bcd.layout.annotations = []
    #     return dcc.Graph(clear_on_unhover=True,
    #                      config={'displayModeBar': False, 'scrollZoom': False},
    #                      id="id_graph7",
    #                      figure=mechanism_dimensions_chart_abc), \
    #         dcc.Graph(clear_on_unhover=True,
    #                   config={'displayModeBar': False, 'scrollZoom': False},
    #                   id="id_graph8",
    #                   figure=mechanism_dimensions_chart_bcd)
    # else:
    #     return (mechanism_dimensions_chart_abc, mechanism_dimensions_chart_bcd)


# def build_mechanism_plot(joint_trajectory, AB, BC, CD):
#     printinfo("build_mechanism_plot", "-")
#     trajectory_path = joint_trajectory
#     yrange = [-(AB + BC + CD), 20]
#     xrange = [-(AB + BC + CD) / 2, (AB + BC + CD) + 20]
#     return dcc.Graph(
#                 config={'displayModeBar': False, 'scrollZoom': False},
#                 id='id_graph3',
#                 figure=go.Figure(
#                     data=[go.Scatter(
#                         x=[joint_trajectory["Ax"][0], joint_trajectory["Bx"][0],
#                            joint_trajectory["Gx"][0]],
#                         y=[joint_trajectory["Ay"][0], joint_trajectory["By"][0],
#                            joint_trajectory["Gy"][0]],
#                         line=dict(color='blue', dash='solid', width=2),
#                         # hovertemplate=None, hoverinfo='skip'
#                     ),
#                         go.Scatter(
#                             x=[joint_trajectory["Fx"][0], joint_trajectory["Bx"][0],
#                                joint_trajectory["Cx"][0]],
#                             y=[joint_trajectory["Fy"][0], joint_trajectory["By"][0],
#                                joint_trajectory["Cy"][0]],
#                             line=dict(color='green', dash='solid', width=2),
#                             # hovertemplate=None, hoverinfo='skip'
#                         ),
#                         go.Scatter(
#                             x=[joint_trajectory["Hx"][0], joint_trajectory["Cx"][0],
#                                joint_trajectory["Dx"][0]],
#                             y=[joint_trajectory["Hy"][0], joint_trajectory["Cy"][0],
#                                joint_trajectory["Dy"][0]],
#                             line=dict(color='orange', dash='solid', width=2),
#                             # hovertemplate=None, hoverinfo='skip'
#                         ),
#                         go.Scatter(
#                             x=[joint_trajectory["Gx"][0], joint_trajectory["Hx"][0]],
#                             y=[joint_trajectory["Gy"][0], joint_trajectory["Hy"][0]],
#                             line=dict(color='purple', dash='solid', width=2),
#                             # hovertemplate=None, hoverinfo='skip'
#                         ),
#                         go.Scatter(
#                             x=[joint_trajectory["Ex"][0], joint_trajectory["Fx"][0]],
#                             y=[joint_trajectory["Ey"][0], joint_trajectory["Fy"][0]],
#                             line=dict(color='purple', dash='solid', width=2),
#                             # hovertemplate=None, hoverinfo='skip'
#                         ),
#
#                         go.Scatter(x=trajectory_path["Ax"], y=trajectory_path["Ay"],
#                                    line=dict(color='gray', dash='dot', width=1),
#                                    hovertemplate=None, hoverinfo='skip'),
#                         go.Scatter(x=trajectory_path["Bx"], y=trajectory_path["By"],
#                                    line=dict(color='gray', dash='dot', width=1),
#                                    hovertemplate=None, hoverinfo='skip'),
#                         go.Scatter(x=trajectory_path["Cx"], y=trajectory_path["Cy"],
#                                    line=dict(color='gray', dash='dot', width=1),
#                                    hovertemplate=None, hoverinfo='skip'),
#                         go.Scatter(x=trajectory_path["Dx"], y=trajectory_path["Dy"],
#                                    line=dict(color='gray', dash='dot', width=1),
#                                    hovertemplate=None, hoverinfo='skip'),
#                         go.Scatter(x=trajectory_path["Ax"], y=trajectory_path["Ay"],
#                                    line=dict(color='gray', dash='dot', width=1),
#                                    hovertemplate=None, hoverinfo='skip'),
#                         go.Scatter(x=trajectory_path["Bx"], y=trajectory_path["By"],
#                                    line=dict(color='gray', dash='dot', width=1),
#                                    hovertemplate=None, hoverinfo='skip'),
#                         go.Scatter(x=trajectory_path["Cx"], y=trajectory_path["Cy"],
#                                    line=dict(color='gray', dash='dot', width=1),
#                                    hovertemplate=None, hoverinfo='skip'),
#                         go.Scatter(x=trajectory_path["Dx"], y=trajectory_path["Dy"],
#                                    line=dict(color='gray', dash='dot', width=1),
#                                    hovertemplate=None, hoverinfo='skip'), ],
#
#                     layout=go.Layout(
#                         hovermode=None,
#                         hoverlabel=None,
#                         # height=350,
#                         # width = 420,
#                         showlegend=False,
#                         margin=dict(
#                             l=0,
#                             r=0,
#                             t=0,
#                             b=0,
#                         ),
#
#                         paper_bgcolor='rgba(0,0,0,0)',
#                         plot_bgcolor='rgba(0,0,0,0)',
#                         # paper_bgcolor='blue',
#                         # plot_bgcolor='pink',
#                         yaxis=dict(
#                             mirror=True,
#                             # tick0 = -80,
#                             # dtick = 20,
#                             fixedrange=True,
#                             range=yrange,
#                             # color = 'none',
#                             # gridcolor = '#103d52',
#                             # gridwidth = 0,
#                             scaleanchor="x",
#                             scaleratio=1,
#                             # linecolor = 'white',
#                             zeroline=False,
#                             color='rgba(0,0,0,0)',
#                             # tickvals = [-80,-60,-40,-20,0,20]
#                             showgrid=False,
#                         ),
#                         xaxis=dict(
#                             showgrid=False,
#                             mirror=True,
#                             # tick0 = -60,
#                             # dtick = 20,
#                             range=xrange,
#                             fixedrange=True,
#
#                             autorange=False,
#                             # gridcolor='#103d52',
#                             # gridwidth=1,
#                             color='rgba(0,0,0,0)',
#                             # linecolor = 'white',
#                             zeroline=False,
#                             scaleanchor="y",
#                             scaleratio=1,
#                         ),
#                         updatemenus=[dict(
#                             # bgcolor = "#103d52",
#                             # bordercolor = '#51656e',
#                             # font = {"color": "#ffd64f"},
#                             direction='left',
#                             # pad = {"r": 10, "t": 87},
#                             showactive=False,
#                             # x = 0.2,
#                             # y = 0.1,
#                             x=0,
#                             y=0.05,
#                             xanchor='left',
#                             yanchor='bottom',
#                             type="buttons",
#                             buttons=[
#                                 dict(
#                                     label="\u23F5",
#                                     method="animate",
#                                     args=[
#                                         None,
#                                         {
#                                             "frame": {
#                                                 "duration": 20,
#                                                 "redraw": True
#                                             },
#                                             "fromcurrent": True,
#                                             "transition": {
#                                                 "duration": 100
#                                             },
#                                             "mode": 'immediate'
#                                         }
#                                     ],
#                                 ),
#                                 dict(
#                                     label="\u23F8",
#                                     method="animate",
#                                     args=[
#                                         [None],
#                                         {
#                                             "frame": {
#                                                 "duration": 0,
#                                                 "redraw": False
#                                             },
#                                             "mode": 'immediate',
#                                             "transition": {
#                                                 "duration": 0
#                                             }
#                                         }
#                                     ],
#                                 )
#                             ]
#                         )],
#
#                         sliders=[dict(steps=[dict(method='animate',
#                                                   args=[[f'frame{k}'],
#                                                         dict(mode='e',
#                                                              frame=dict(duration=0, redraw=False),
#                                                              transition=dict(duration=100))
#                                                         ],
#                                                   # label=f'\u03B1 = {trajectory_path["alpha"][k]}\u00B0'
#                                                   # label = f'{k/4}',
#                                                   label="{:.2f}".format(k / 4),
#                                                   # label = ' ',
#                                                   # label = [1,2,3],
#
#                                                   visible=True,
#                                                   # label = f'\u03B1 = {trajectory_path["alpha"][k]}\u00B0 tet1 = {round(trajectory_path["tet1"][k],2)} gamCF = {round(trajectory_path["gamCF"][k],2)}'
#                                                   ) for k in range(len(trajectory_path["alpha"]))],
#                                       transition=dict(duration=100, easing="cubic-in-out"),
#                                       x=0.25,
#                                       y=0,
#                                       # pad = {"b": 10, "t": 50},
#                                       yanchor="bottom",
#                                       # xanchor = "left",
#                                       # tickcolor = 'red',
#                                       font={"color": "rgba(0,0,0,0)"},
#                                       # direction = 'up',
#                                       currentvalue=dict(font=dict(size=12, color="#ffd64f"),
#                                                         prefix='\u03B1 = ',
#                                                         suffix=' \u00B0',
#                                                         visible=True,
#                                                         xanchor='right',
#                                                         # wyrównanie aktualnej wartości do prawej strony slidera
#                                                         # yanchor='top'
#
#                                                         ),
#                                       len=0.6,
#                                       active=0,
#                                       minorticklen=0,
#                                       ticklen=0,
#                                       borderwidth=5,
#                                       bordercolor='green',
#                                       bgcolor='pink',
#                                       activebgcolor='blue',
#
#                                       )
#                                  ]
#                     ),
#                     frames=[
#                         go.Frame(
#                             data=[
#                                 go.Scatter(
#                                     x=[joint_trajectory["Ax"][k], joint_trajectory["Bx"][k],
#                                        joint_trajectory["Gx"][k]],
#                                     y=[joint_trajectory["Ay"][k], joint_trajectory["By"][k],
#                                        joint_trajectory["Gy"][k]],
#                                     line=dict(color='blue', dash='solid', width=4),
#                                     hovertemplate=None
#                                 ),
#                                 go.Scatter(
#                                     x=[joint_trajectory["Fx"][k], joint_trajectory["Bx"][k],
#                                        joint_trajectory["Cx"][k]],
#                                     y=[joint_trajectory["Fy"][k], joint_trajectory["By"][k],
#                                        joint_trajectory["Cy"][k]],
#                                     line=dict(color='green', dash='solid', width=4),
#                                     hovertemplate=None
#                                 ),
#                                 go.Scatter(
#                                     x=[joint_trajectory["Hx"][k], joint_trajectory["Cx"][k],
#                                        joint_trajectory["Dx"][k]],
#                                     y=[joint_trajectory["Hy"][k], joint_trajectory["Cy"][k],
#                                        joint_trajectory["Dy"][k]],
#                                     line=dict(color='orange', dash='solid', width=4),
#                                     hovertemplate=None
#                                 ),
#                                 go.Scatter(
#                                     x=[joint_trajectory["Gx"][k], joint_trajectory["Hx"][k]],
#                                     y=[joint_trajectory["Gy"][k], joint_trajectory["Hy"][k]],
#                                     line=dict(color='purple', dash='solid', width=2),
#                                     hovertemplate=None
#                                 ),
#                                 go.Scatter(
#                                     x=[joint_trajectory["Ex"][k], joint_trajectory["Fx"][k]],
#                                     y=[joint_trajectory["Ey"][k], joint_trajectory["Fy"][k]],
#                                     line=dict(color='purple', dash='solid', width=2),
#                                     hovertemplate=None
#                                 )
#                             ],
#                             name=f'frame{k}'
#                         )
#                         for k in range(len(joint_trajectory["Ex"]))]
#
#                 )
#             )


def build_mechanism_position_chart(type, joint_trajectory, ab, bc, cd):
    # The function builds graphs "mechanism-positions-plot" and "mechanism-vectors-plot".
    # The function is called two times by "manage_content_style_in_tab1()"

    # Init basics variables
    yrange = [-(ab + bc + cd), 20]
    xrange = [-(ab + bc + cd) / 2, (ab + bc + cd) + 20]

    # Prepare the keys used to create the list of traces (for x and y in "prepare_position_traces()").
    # "key_names_for_mechanism_lines" is necessary to build both graphs.
    key_names_for_mechanism_lines = [
        [["Ax", "Bx", "Gx"], ["Ay", "By", "Gy"]], [["Fx", "Bx", "Cx"], ["Fy", "By", "Cy"]],
        [["Hx", "Cx", "Dx"], ["Hy", "Cy", "Dy"]], [["Gx", "Hx"], ["Gy", "Hy"]],
        [["Ex", "Fx"], ["Ey", "Fy"]], ["Bx", "By"], ["Cx", "Cy"], ["Dx", "Dy"]]

    def prepare_position_traces():
        traces_pos = []
        for trace_no, key_names in enumerate(key_names_for_mechanism_lines):
            if trace_no < 5:
                x = [joint_trajectory[element][0] for element in key_names[0]]
                y = [joint_trajectory[element][0] for element in key_names[1]]
            else:
                x = joint_trajectory[key_names[0]]
                y = joint_trajectory[key_names[1]]
            traces_pos.append(
                go.Scatter(
                    x=x, y=y, line=mechanism_lines[trace_no],
                    marker=dict(opacity=1, line=dict(width=1, color='#383734')),
                    hovertemplate=hovertemplates[trace_no], hoverinfo=hoverinfos[trace_no]))
        return traces_pos

    def prepare_vectors_traces():
        traces_vec = []
        for trace_no, key_names in enumerate(key_names_for_vector_lines):
            for i in range(2):
                traces_vec.append(
                    go.Scatter(
                        x=[joint_trajectory[element][0] for element in key_names[i][0]],
                        y=[joint_trajectory[element][0] for element in key_names[i][1]],
                        line=vector_lines[trace_no][i], marker=dict(
                            symbol=joint_trajectory[key_names_for_vector_markers[trace_no][i]][0],
                            color=vector_marker_colors[trace_no], line=dict(color="black", width=1), size=[3, 6]),
                        hovertemplate=hovertemplates[trace_no], hoverinfo=hoverinfos[trace_no]))
        return traces_vec

    def set_layout():
        return go.Layout(hovermode=None, hoverlabel=None, showlegend=False,
                         margin=dict(l=0, r=0, t=0, b=0, pad=0),
                         paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)', autosize=True,
                         yaxis=dict(fixedrange=True, range=yrange, scaleanchor="x", scaleratio=1,
                                    zeroline=False, color='rgba(0,0,0,0)', showgrid=False, visible=False),
                         xaxis=dict(showgrid=False, range=xrange, fixedrange=True, autorange=False, visible=False,
                                    color='rgba(0,0,0,0)', zeroline=False, scaleanchor="y", scaleratio=1))

    if type == "normal":
        chart_id_name = "mechanism-positions-plot"

        # Prepare the lines of the mechanism.
        mechanism_lines = [
            dict(color='blue', dash='solid', width=2), dict(color='green', dash='solid', width=2),
            dict(color='orange', dash='solid', width=2), dict(color='purple', dash='solid', width=2),
            dict(color='purple', dash='solid', width=2), dict(color='gray', dash='dot', width=1),
            dict(color='gray', dash='dot', width=1), dict(color='gray', dash='dot', width=1)]

        # TODO - prepare hovertemplates
        hovertemplates = [None] * 8
        hoverinfos = ["skip"] * 8

        traces = prepare_position_traces()

    elif type == "vectors":
        chart_id_name = "mechanism-vectors-plot"

        # Prepare the lines of the mechanism.
        mechanism_line_dash = ['solid', 'solid', 'solid', 'solid', 'solid', 'dot', 'dot', 'dot']
        mechanism_lines = [dict(color='gray', dash=dash, width=1) for dash in mechanism_line_dash]

        # Prepare vectors line style. "vector_lines" is a 6-element list (6 points), where each element contain two
        # dictionaries describing the vector lines in the x and y direction
        vector_line_colors = ['rgba(240,98,108,0.9)', 'rgba(240,98,108,0.9)', 'rgba(240,98,108,0.9)',
                              'rgba(240, 216, 84,0.9)', 'rgba(240, 216, 84,0.9)', 'rgba(240, 216, 84,0.9)']
        vector_lines = [[dict(color=color, dash='solid', width=1)] * 2 for color in vector_line_colors]

        # Prepare the colors of the vector markers. There is a list with six elements. Each element is list has two
        # values of color - one for points B, C, D, F, G, H and one for end of vectors.
        vector_marker_colors = \
            [['rgba(0,0,0,0)', 'rgba(240,98,108,0.9)']] * 3 + [['rgba(0,0,0,0)', 'rgba(240, 216, 84,0.9)']] * 3

        # TODO - prepare hovertemplates
        hovertemplates = [None] * 12
        hoverinfos = ["skip"] * 12

        # Prepare the keys used to create the list of traces (for x and y in "prepare_vectors_traces()").
        key_names_for_vector_lines = [
            [[["Bx", "vhBx"], ["By", "vhBy"]], [["Bx", "vvBx"], ["By", "vvBy"]]],
            [[["Cx", "vhCx"], ["Cy", "vhCy"]], [["Cx", "vvCx"], ["Cy", "vvCy"]]],
            [[["Dx", "vhDx"], ["Dy", "vhDy"]], [["Dx", "vvDx"], ["Dy", "vvDy"]]],
            [[["Fx", "vhFx"], ["Fy", "vhFy"]], [["Fx", "vvFx"], ["Fy", "vvFy"]]],
            [[["Gx", "vhGx"], ["Gy", "vhGy"]], [["Gx", "vvGx"], ["Gy", "vvGy"]]],
            [[["Hx", "vhHx"], ["Hy", "vhHy"]], [["Hx", "vvHx"], ["Hy", "vvHy"]]]]

        # Prepare the keys used to create the list of traces (for markers color in "prepare_vectors_traces()").
        key_names_for_vector_markers = [
            ["vhsB", "vvsB"], ["vhsC", "vvsC"], ["vhsD", "vvsD"], ["vhsF", "vvsF"], ["vhsG", "vvsG"], ["vhsH", "vvsH"]]

        trace_pos = prepare_position_traces()
        trace_vek = prepare_vectors_traces()
        traces = trace_pos + trace_vek

    return dcc.Graph(
        config={'displayModeBar': False, 'scrollZoom': False}, id=chart_id_name,
        figure=go.Figure(data=traces, layout=set_layout()))


def build_vector_velocity_chart(joint_trajectory, AB, BC, CD):
    # Function returned velocity vectors in one of charts in tab "model parameters"
    # There are many traces, so it is executing in two for loops based on the prepared dict with parameters
    printinfo("build_mechanism_plot", "-")
    trajectory_path = joint_trajectory
    yrange = [-(AB + BC + CD), 20]
    xrange = [-(AB + BC + CD) / 2, (AB + BC + CD) + 20]

    def set_data(joint_trajectory):
        data = []
        parameters_dict_mechanism = dict(
            trace1=dict(x=[joint_trajectory["Ax"][0], joint_trajectory["Bx"][0], joint_trajectory["Gx"][0]],
                        y=[joint_trajectory["Ay"][0], joint_trajectory["By"][0], joint_trajectory["Gy"][0]],
                        dash='solid'),
            trace2=dict(x=[joint_trajectory["Fx"][0], joint_trajectory["Bx"][0], joint_trajectory["Cx"][0]],
                        y=[joint_trajectory["Fy"][0], joint_trajectory["By"][0], joint_trajectory["Cy"][0]],
                        dash='solid'),
            trace3=dict(x=[joint_trajectory["Hx"][0], joint_trajectory["Cx"][0], joint_trajectory["Dx"][0]],
                        y=[joint_trajectory["Hy"][0], joint_trajectory["Cy"][0], joint_trajectory["Dy"][0]],
                        dash='solid'),
            trace4=dict(x=[joint_trajectory["Gx"][0], joint_trajectory["Hx"][0]],
                        y=[joint_trajectory["Gy"][0], joint_trajectory["Hy"][0]], dash='solid'),
            trace5=dict(x=[joint_trajectory["Ex"][0], joint_trajectory["Fx"][0]],
                        y=[joint_trajectory["Ey"][0], joint_trajectory["Fy"][0]], dash='solid'),
            trace6=dict(x=trajectory_path["Bx"], y=trajectory_path["By"], dash='dot'),
            trace7=dict(x=trajectory_path["Cx"], y=trajectory_path["Cy"], dash='dot'),
            trace8=dict(x=trajectory_path["Dx"], y=trajectory_path["Dy"], dash='dot'))

        for trace in parameters_dict_mechanism.values():
            data.append(go.Scatter(x=trace["x"], y=trace["y"], line=dict(color='gray', dash=trace["dash"], width=1)))

        parameters_dict_vectors = dict(
            trace1=dict(x=[trajectory_path["Bx"][0], trajectory_path["Bx"][0] + trajectory_path["dBx"][0]],
                        y=[trajectory_path["By"][0], trajectory_path["By"][0]],
                        symbol=["circle", "diamond-wide"], color='rgba(240,98,108,0.5)'),
            trace2=dict(x=[trajectory_path["Cx"][0], trajectory_path["Cx"][0] + trajectory_path["dCx"][0]],
                        y=[trajectory_path["Cy"][0], trajectory_path["Cy"][0]],
                        symbol=["circle", "diamond-wide"], color='rgba(240,98,108,0.5)'),
            trace3=dict(x=[trajectory_path["Dx"][0], trajectory_path["Dx"][0] + trajectory_path["dDx"][0]],
                        y=[trajectory_path["Dy"][0], trajectory_path["Dy"][0]],
                        symbol=["circle", "diamond-wide"], color='rgba(240,98,108,0.5)'),
            trace4=dict(x=[trajectory_path["Fx"][0], trajectory_path["Fx"][0] + trajectory_path["dFx"][0]],
                        y=[trajectory_path["Fy"][0], trajectory_path["Fy"][0]],
                        symbol=["circle", "diamond-wide"], color='rgba(121,175,232,0.5)'),
            trace5=dict(x=[trajectory_path["Gx"][0], trajectory_path["Gx"][0] + trajectory_path["dGx"][0]],
                        y=[trajectory_path["Gy"][0], trajectory_path["Gy"][0]],
                        symbol=["circle", "diamond-wide"], color='rgba(121,175,232,0.5)'),
            trace6=dict(x=[trajectory_path["Hx"][0], trajectory_path["Hx"][0] + trajectory_path["dHx"][0]],
                        y=[trajectory_path["Hy"][0], trajectory_path["Hy"][0]],
                        symbol=["circle", "diamond-wide"], color='rgba(121,175,232,0.5)'),
            trace7=dict(x=[trajectory_path["Bx"][0], trajectory_path["Bx"][0]],
                        y=[trajectory_path["By"][0], trajectory_path["By"][0] + trajectory_path["dBy"][0]],
                        symbol=["circle", "diamond-tall"], color='rgba(240,98,108,0.5)'),
            trace8=dict(x=[trajectory_path["Cx"][0], trajectory_path["Cx"][0]],
                        y=[trajectory_path["Cy"][0], trajectory_path["Cy"][0] + trajectory_path["dCy"][0]],
                        symbol=["circle", "diamond-tall"], color='rgba(240,98,108,0.5)'),
            trace9=dict(x=[trajectory_path["Dx"][0], trajectory_path["Dx"][0]],
                        y=[trajectory_path["Dy"][0], trajectory_path["Dy"][0] + trajectory_path["dDy"][0]],
                        symbol=["circle", "diamond-tall"], color='rgba(240,98,108,0.5)'),
            trace10=dict(x=[trajectory_path["Fx"][0], trajectory_path["Fx"][0]],
                         y=[trajectory_path["Fy"][0], trajectory_path["Fy"][0] + trajectory_path["dFy"][0]],
                         symbol=["circle", "diamond-tall"], color='rgba(121,175,232,0.5)'),
            trace11=dict(x=[trajectory_path["Gx"][0], trajectory_path["Gx"][0]],
                         y=[trajectory_path["Gy"][0], trajectory_path["Gy"][0] + trajectory_path["dGy"][0]],
                         symbol=["circle", "diamond-tall"], color='rgba(121,175,232,0.5)'),
            trace12=dict(x=[trajectory_path["Hx"][0], trajectory_path["Hx"][0]],
                         y=[trajectory_path["Hy"][0], trajectory_path["Hy"][0] + trajectory_path["dHy"][0]],
                         symbol=["circle", "diamond-tall"], color='rgba(121,175,232,0.5)'))

        for trace in parameters_dict_vectors.values():
            data.append(go.Scatter(x=trace["x"], y=trace["y"], line=dict(color=trace["color"], dash='solid', width=1),
                                   marker=dict(opacity=0.7, color=trace["color"], line=dict(color='white', width=0),
                                               symbol=trace["symbol"], size=[4, 6])))
        return data

    def set_layout():
        return go.Layout(hovermode=None, hoverlabel=None, showlegend=False, margin=dict(l=0, r=0, t=0, b=0),
                         paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)',
                         yaxis=dict(mirror=True, fixedrange=True, range=yrange, scaleanchor="x", scaleratio=1,
                                    zeroline=False, color='rgba(0,0,0,0)', showgrid=False),
                         xaxis=dict(showgrid=False, mirror=True, range=xrange, fixedrange=True, autorange=False,
                                    color='rgba(0,0,0,0)', zeroline=False, scaleanchor="y", scaleratio=1))

    return dcc.Graph(
        config={'displayModeBar': False, 'scrollZoom': False}, id='gfgf',
        figure=go.Figure(data=set_data(joint_trajectory), layout=set_layout()))


def build_bar_chart(direction, joint_trajectory):
    # The function builds bar chart. First, the data is prepared ("sufix", lists, "color" etc). Next, chart is building
    # based on "set_data()" and "set_layout()"
    point_names = ["B", "C", "D", "F", "G", "H"]
    if direction == "x":
        sufix, id = "x", 'id_graph5'
        parameters = [point_name + "_x" for point_name in point_names]
        color, title = "#6be0e8", r"$x\text{-direction velocity  } {[\frac{\text{mm}}{^\circ}]}$"
    elif direction == "y":
        sufix, id = "y", 'id_graph6'
        parameters = [point_name + "_y" for point_name in point_names]
        color, title = "#f7d479", r"$y\text{-direction velocity  } {[\frac{\text{mm}}{^\circ}]}$"
    elif direction == "resultant":
        sufix, id = "", 'id_graph11'
        parameters = point_names
        color, title = "#acf279", r"$\text{resultant velocity  } {[\frac{\text{mm}}{^\circ}]}$"

    def set_data(joint_trajectory):
        labels = [r"${\frac{\Delta " + param + r"}{\Delta\alpha}}$" for param in parameters]
        values = [joint_trajectory["d" + point_name + sufix][0] for point_name in point_names]
        return go.Bar(x=labels, y=np.abs(values), text=values, textposition="outside",
                      outsidetextfont=dict(color=color), marker=dict(color=color), hovertemplate=None, hoverinfo='skip')

    def set_layout():
        return go.Layout(
            hovermode=None, hoverlabel=None, showlegend=False, margin=dict(l=10, r=10, t=20, b=10),
            paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)',
            title=dict(text=title, font=dict(color="white", size=12), xanchor='center', x=0.5, xref='paper',
                       yanchor='bottom',
                       y=0.95, yref='paper'),
            yaxis=dict(mirror=True, fixedrange=True, showticklabels=False, showgrid=False, range=[-0.05, 1.4],
                       zeroline=False, color='white'),
            xaxis=dict(fixedrange=True, showgrid=False, mirror=True, ticklabelposition="outside top", color='white',
                       zeroline=True))

    return dcc.Graph(config={'displayModeBar': False, 'scrollZoom': False}, id=id, mathjax=True,
                     figure=go.Figure(data=set_data(joint_trajectory), layout=set_layout()))


def build_slider_mechanism_params_content(jt_as_dict):
    step = 0.1  # check update_slider_value_model_parameters, update_model_parameters_plot calculate_model
    return [dcc.Slider(0, 90, step, id="slider", updatemode='drag', value=0, marks=None,
                       tooltip={"placement": "bottom", "always_visible": True}), \
            html.Div(id="button_container", children=[
                # label="\u23F8", label = "\u23F5",
                dbc.Button("\u23F5", id="button_model_parameters_start",
                           className="my-button-model-parameters", n_clicks=0),
                dbc.Button("\u23F8", id="button_model_parameters_pause",
                           className="my-button-model-parameters", n_clicks=0),
                dbc.Button("\u23ED", id="button_model_parameters_next",
                           className="my-button-model-parameters", n_clicks=0)
            ])]


def build_static_velocity_chart(joint_trajectory):
    max_value = np.max([np.max(joint_trajectory[key]) for key in ["dB", "dC", "dD", "dF", "dG", "dH"]])
    y_range = [0, 1.3 * max_value]
    title = r"$\text{  resultant velocity}$"
    yticks = np.round(np.linspace(0, max_value * 1.3, 5), 3)
    xticks = [0, 15, 30, 45, 60, 75, 90, 100]

    figure = go.Figure(
        data=[go.Scatter(x=joint_trajectory["alpha"][1:], y=joint_trajectory[name][1:],
                         name=r"$" + name[1].lower() + r"(\alpha)$", hovertemplate=None, hoverinfo="skip",
                         mode='lines', marker=dict(symbol="circle", size=1, line=dict(width=0)))
              for idx, name in enumerate(["dB", "dC", "dD", "dF", "dG", "dH"])],
        layout=dict(
            hovermode=None,
            title=dict(text=title, font=dict(color="white", size=12), xanchor='center', x=0.5, xref='paper',
                       yanchor='bottom',
                       y=1.0, yref='paper'),
            legend=dict(xanchor="left", yanchor="middle", x=1, y=0.5, orientation="v",
                        font=dict(color='white', size=9), itemwidth=30),
            margin=dict(l=10, r=10, t=20, b=10), paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)',
            xaxis=dict(color="white", zeroline=False, range=[0, 100], showgrid=False, automargin=True, mirror=False,
                       tickfont=dict(size=10), ticks='inside', showline=True, linecolor="gray", linewidth=1,
                       ticklabelposition="outside", tickcolor="#ffd64f", tickvals=xticks,
                       ticktext=[f"${tick_value:.0f}$" if idx < 7
                                 else r"$[^\circ]$"
                                 for idx, tick_value in enumerate(xticks)]),
            yaxis=dict(color="white", zeroline=False, range=y_range, showgrid=False, griddash='solid', gridcolor='gray',
                       gridwidth=1, tickfont=dict(size=10, color="white"), mirror=False, ticks='inside', showline=True,
                       tickvals=yticks, ticklabelposition="outside", linecolor="gray", tickcolor="#ffd64f",
                       ticktext=[r"$" + f"{tick_value:.3f}$" if idx < 4
                                 else r"$[\frac{\text{mm}}{^\circ}]$"
                                 for idx, tick_value in enumerate(yticks)])))
    return [
        html.Div(id='id21', children=[
            dcc.Graph(id="mechanism-velocities-static", config={'displayModeBar': False, 'scrollZoom': False},
                      mathjax=True, figure=figure)])]


def build_static_angles_chart(joint_trajectory):
    max_value = np.max([np.max(joint_trajectory[key]) for key in ["dB", "dC", "dD", "dF", "dG", "dH"]])
    title = r"$$\text{angles between segments}$$"

    yticks = [90, 105, 120, 135, 150, 165, 180, 210]
    xticks = [0, 15, 30, 45, 60, 75, 90, 100]

    trace_names = [r"$\gamma_1$", r"$\gamma_2$"]
    figure = go.Figure(
        data=[go.Scatter(x=joint_trajectory["alpha"], y=joint_trajectory[key],
                         name=trace_names[idx],
                         hovertemplate=None, hoverinfo="skip",
                         mode='lines')
              for idx, key in enumerate(["kABC", "kBCD"])],
        layout=dict(
            # hovermode='x unified',
            title=dict(text=title, font=dict(color="white", size=12), xanchor='center', x=0.5, xref='paper',
                       yanchor='bottom',
                       y=1, yref='paper'),
            legend=dict(xanchor="right", yanchor="middle", x=1, y=0.5, orientation="v",
                        font=dict(color='white', size=9), itemwidth=30),
            margin=dict(l=40, r=0, t=20, b=20), paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)',
            xaxis=dict(color="white", zeroline=False, showgrid=False, automargin=True, mirror=False, gridcolor='gray',
                       gridwidth=1, griddash='dot', range=[0, 100],
                       tickfont=dict(size=10), ticks='inside', showline=True, linecolor="gray", linewidth=1,
                       ticklabelposition="outside", tickcolor="#ffd64f",
                       tickvals=xticks,
                       ticktext=[f"${tick_value:.0f}$" if idx < 7
                                 else r"$[^\circ]$"
                                 for idx, tick_value in enumerate(xticks)]
                       ),
            yaxis=dict(color="white", zeroline=False, showgrid=False, griddash='dot', gridcolor='gray',
                       gridwidth=1, tickfont=dict(size=10, color="white"), mirror=False, ticks='inside', showline=True,
                       ticklabelposition="outside", linecolor="gray", range=[70, 210], tickcolor="#ffd64f",
                       tickvals=yticks,
                       ticktext=[f"${tick_value:.0f}$" if idx < 7
                                 else r"$[^\circ]$"
                                 for idx, tick_value in enumerate(yticks)]
                       )))
    #
    color = ["rgb(99, 110, 250)", "rgb(239, 85, 59)"]
    text_color = ["white", "white"]
    values = [round(joint_trajectory[name][0], 2) for name in ["kABC", "kBCD"]]
    text_values = [f"{value:.2f}\u00B0" for value in values]
    bar_figure = go.Figure(
        data=go.Bar(x=trace_names, y=90 + np.abs(values), text=text_values, textposition="inside", textangle=-90,
                    textfont=dict(size=13, color="white"), hovertemplate=None, hoverinfo="skip",
                    marker=dict(color=color)),
        layout=go.Layout(
            hovermode=None, hoverlabel=None, showlegend=False, margin=dict(l=0, r=10, t=20, b=20),
            paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)',
            title=dict(text=r"$\alpha\text{ = } 180^\circ$", font=dict(color="white", size=12),
                       xanchor='center', x=0.5, xref='paper', yanchor='top',
                       y=1, yref='paper'),
            yaxis=dict(mirror=True, fixedrange=True, showticklabels=False, showgrid=False, range=[0, 300],
                       zeroline=False, color='white'),
            xaxis=dict(fixedrange=True, showgrid=False, mirror=True, ticklabelposition="outside top", color='white',
                       zeroline=True))
    )

    return [
        html.Div(id='id22', children=[
            dcc.Graph(id="mechanism-angles-static", config={'displayModeBar': False, 'scrollZoom': False},
                      mathjax=True, figure=figure)]),
        html.Div(id='id23', children=[
            dcc.Graph(id="mechanism-angles-dynamic", config={'displayModeBar': False, 'scrollZoom': False},
                      mathjax=True, figure=bar_figure)])
    ]


def build_description_chart():

    ang1, ang2, ang3 = 20, 30, 45  # the inclination angles of the segments for alpha
    ang1_prim, ang2_prim, ang3_prim = 30, 60, 75  # the inclination angles of the segments for alpha prim
    ang1_aend, ang2_aend, ang3_aend = 85, 170, 260  # the inclination angles of the segments for almost end alpha

    # position of points for alpha
    ax, ay = 0, 0
    bx, by = ax + 40 * np.cos(ang1 * np.pi / 180), ay - 40 * np.sin(ang1 * np.pi / 180)
    cx, cy = bx + 35 * np.cos(ang2 * np.pi / 180), by - 35 * np.sin(ang2 * np.pi / 180)
    dx, dy = cx + 25 * np.cos(ang3 * np.pi / 180), cy - 25 * np.sin(ang3 * np.pi / 180)
    # position of points for alpha prim
    bx_prim, by_prim = ax + 40 * np.cos(ang1_prim * np.pi / 180), ay - 40 * np.sin(ang1_prim * np.pi / 180)
    cx_prim, cy_prim = bx_prim + 35 * np.cos(ang2_prim * np.pi / 180), by_prim - 35 * np.sin(ang2_prim * np.pi / 180)
    dx_prim, dy_prim = cx_prim + 25 * np.cos(ang3_prim * np.pi / 180), cy_prim - 25 * np.sin(ang3_prim * np.pi / 180)
    # position of points for almost end alpha
    bx_aend, by_aend = ax + 40 * np.cos(ang1_aend * np.pi / 180), ay - 40 * np.sin(ang1_aend * np.pi / 180)
    cx_aend, cy_aend = bx_aend + 35 * np.cos(ang2_aend * np.pi / 180), by_aend - 35 * np.sin(ang2_aend * np.pi / 180)
    dx_aend, dy_aend = cx_aend + 25 * np.cos(ang3_aend * np.pi / 180), cy_aend - 25 * np.sin(ang3_aend * np.pi / 180)

    r1, r2, r3, r4 = 70, 50, 60, 30  # length of extension lines r1: alpha zero, r2: alpha, r3: alpha prim, r4: gammas
    arc_offset = 15  # additional arc span in degree
    radius_offset = 5  # dimensional arc offset in millimeters
    arc_offset_degree_base = 15  # 15 degree for arc radius equal 70 mm
    arc_offset_radius_base = 70  # 15 degree for arc radius equal 70 mm
    # calculate arc offset for each case // arc offset - additional arc span in degrees
    arc_offset_alpha = \
        2 * np.arcsin(arc_offset_radius_base * np.sin(arc_offset_degree_base / 2 * np.pi / 180) /
                      (r2 - radius_offset)) * 180 / np.pi
    arc_offset_alpha_prim = \
        2 * np.arcsin(arc_offset_radius_base * np.sin(arc_offset_degree_base / 2 * np.pi / 180) /
                      (r3 - radius_offset)) * 180 / np.pi
    arc_offset_alpha_aend = \
        2 * np.arcsin(arc_offset_radius_base * np.sin(arc_offset_degree_base / 2 * np.pi / 180) /
                      (r4 - radius_offset)) * 180 / np.pi
    # convert to int (for loop range required integer)
    arc_offset_alpha = int(round(arc_offset_alpha, 0))
    arc_offset_alpha_prim = int(round(arc_offset_alpha_prim, 0))

    # dictionary list containing coordinates of dimensional arcs
    # curves_alpha for alpha, curves_alpha_prim for alpha prim
    # dimensional arc divided into two fragments -> two elements in the list
    curves_alpha = [
        dict(x=[ax - (r2 - radius_offset) * np.cos(beta * np.pi / 180) for beta in range(arc_offset_alpha + 1)],
             y=[ay - (r2 - radius_offset) * np.sin(beta * np.pi / 180) for beta in range(arc_offset_alpha + 1)],
             marker=dict(size=[8] + [0 for _beta in range(arc_offset_alpha)], symbol="arrow", angle=0,
                         line=dict(width=0), opacity=1)),
        dict(x=[ax - (r2 - radius_offset) * np.cos(beta * np.pi / 180) for beta in range(ang1, ang1 + arc_offset_alpha + 1)],
             y=[ay + (r2 - radius_offset) * np.sin(beta * np.pi / 180) for beta in range(ang1, ang1 + arc_offset_alpha + 1)],
             marker=dict(size=[8] + [0 for _beta in range(arc_offset_alpha)], symbol="arrow", angle=180 + ang1,
                         line=dict(width=0), opacity=1))]
    curves_alpha_prim = [
        dict(x=[ax - (r3 - radius_offset) * np.cos(beta * np.pi / 180) for beta in range(arc_offset_alpha_prim + 1)],
             y=[ay - (r3 - radius_offset) * np.sin(beta * np.pi / 180) for beta in range(arc_offset_alpha_prim + 1)],
             marker=dict(size=[8] + [0 for _beta in range(arc_offset_alpha_prim)], symbol="arrow", angle=0,
                         line=dict(width=0), opacity=1)),
        dict(x=[ax - (r3 - radius_offset) * np.cos(beta * np.pi / 180) for beta in
                range(ang1_prim, ang1_prim + arc_offset_alpha_prim + 1)],
             y=[ay + (r3 - radius_offset) * np.sin(beta * np.pi / 180) for beta in
                range(ang1_prim, ang1_prim + arc_offset_alpha_prim + 1)],
             marker=dict(size=[8] + [0 for _beta in range(arc_offset_alpha_prim)], symbol="arrow", angle=180 + ang1_prim,
                         line=dict(width=0), opacity=1))]

    # calculate delta D dimension line markers angel
    deltaDx = ((dx - dx_prim) ** 2) ** 0.5
    deltaD = ((dx - dx_prim) ** 2 + (dy - dy_prim) ** 2) ** 0.5
    ws_dim_arrow_angle = np.arcsin(deltaDx / deltaD) * 180 / np.pi  # west-south arrow angle // arcsin(deltaDx/deltaD)

    arrow_offset = 2  # space between point and end of arrow marker
    ext_color = "gray"  # extension line color
    mech_color = "rgba(240,98,108,0.9)"  # mechanism line color
    dim_color = "#ffd64f"  # dimension value color

    extension_lines = [
        dict(x=[ax, ax - r1], y=[ay, ay], code="ext_line01"),
        dict(x=[ax, ax - r2 * np.cos(ang1 * np.pi / 180)],
             y=[ay, ay + r2 * np.sin(ang1 * np.pi / 180)], code="ext_line02"),
        dict(x=[ax, ax - r3 * np.cos(ang1_prim * np.pi / 180)],
             y=[ay, ay + r3 * np.sin(ang1_prim * np.pi / 180)], code="ext_line03"),
        dict(x=[bx_aend, bx_aend + r4 * np.cos(ang1_aend * np.pi / 180)],
             y=[by_aend, by_aend - r4 * np.sin(ang1_aend * np.pi / 180)], code="ext_line04"),
        dict(x=[bx_aend, bx_aend - r4 * np.cos(ang2_aend * np.pi / 180)],
             y=[by_aend, by_aend + r4 * np.sin(ang2_aend * np.pi / 180)], code="ext_line05"),
        dict(x=[cx_aend, cx_aend + r4 * np.cos(ang2_aend * np.pi / 180)],
             y=[cy_aend, cy_aend - r4 * np.sin(ang2_aend * np.pi / 180)], code="ext_line06"),
        dict(x=[cx_aend, cx_aend - r4 * np.cos(ang3_aend * np.pi / 180)],
             y=[cy_aend, cy_aend + r4 * np.sin(ang3_aend * np.pi / 180)], code="ext_line07"),
        dict(x=[dx_prim, dx_prim], y=[dy_prim, dy_prim - 20], code="ext_line09"),
        dict(x=[dx, dx], y=[dy, dy_prim - 20], code="ext_line08"),
        dict(x=[dx_prim, dx + 20], y=[dy_prim, dy_prim], code="ext_line11"),
        dict(x=[dx, dx + 20], y=[dy, dy], code="ext_line10")]

    dimensional_arcs = [
        dict(x=[bx_aend + (r4 - radius_offset) * np.sin(beta * np.pi / 180) for beta in range(5, 100 + 1)],
             y=[by_aend - (r4 - radius_offset) * np.cos(beta * np.pi / 180) for beta in range(5, 100 + 1)],
             code="dim_arc01"),
        dict(x=[cx_aend - (r4 - radius_offset) * np.cos(beta * np.pi / 180) for beta in range(10, 100 + 1)],
             y=[cy_aend - (r4 - radius_offset) * np.sin(beta * np.pi / 180) for beta in range(10, 100 + 1)],
             code="dim_arc02"),
        dict(x=[ax - (r2 - radius_offset) * np.cos(beta * np.pi / 180) for beta in range(arc_offset_alpha + 1)],
             y=[ay - (r2 - radius_offset) * np.sin(beta * np.pi / 180) for beta in range(arc_offset_alpha + 1)],
             code="dim_arc03.1"),
        dict(x=[ax - (r2 - radius_offset) * np.cos(beta * np.pi / 180)
                for beta in range(ang1, ang1 + arc_offset_alpha + 1)],
             y=[ay + (r2 - radius_offset) * np.sin(beta * np.pi / 180)
                for beta in range(ang1, ang1 + arc_offset_alpha + 1)],
             code="dim_arc03.2"),
        dict(x=[ax - (r3 - radius_offset) * np.cos(beta * np.pi / 180) for beta in range(arc_offset_alpha_prim + 1)],
             y=[ay - (r3 - radius_offset) * np.sin(beta * np.pi / 180) for beta in range(arc_offset_alpha_prim + 1)],
             code="dim_arc04.1"),
        dict(x=[ax - (r3 - radius_offset) * np.cos(beta * np.pi / 180)
                for beta in range(ang1_prim, ang1_prim + arc_offset_alpha_prim + 1)],
             y=[ay + (r3 - radius_offset) * np.sin(beta * np.pi / 180)
                for beta in range(ang1_prim, ang1_prim + arc_offset_alpha_prim + 1)],
             code="dim_arc04.2")]

    dimensional_lines = [
        dict(x=[dx - arrow_offset, dx_prim + arrow_offset], y=[dy - arrow_offset, dy_prim + arrow_offset],
             code="dim_line01"),
        dict(x=[dx, dx_prim], y=[dy_prim - 20, dy_prim - 20], code="dim_line02"),
        dict(x=[dx + 20, dx + 20], y=[dy, dy_prim], code="dim_line03")]

    data = [
        # extension line alpha zero X
        go.Scatter(x=[ax, ax - r1], y=[ay, ay], hoverinfo="skip",
                   line=dict(width=1, dash="dot", color=ext_color), mode="lines"),
        # extension line alpha X
        go.Scatter(x=[ax, ax - r2 * np.cos(ang1 * np.pi / 180)], y=[ay, ay + r2 * np.sin(ang1 * np.pi / 180)],
                   hoverinfo="skip", line=dict(width=1, dash="dot", color=ext_color), mode="lines"),
        # extension line alpha prim X
        go.Scatter(x=[ax, ax - r3 * np.cos(ang1_prim * np.pi / 180)], y=[ay, ay + r3 * np.sin(ang1_prim * np.pi / 180)],
                   hoverinfo="skip", line=dict(width=1, dash="dot", color=ext_color), mode="lines"),
        # extension line gamma1 AB X
        go.Scatter(x=[bx_aend, bx_aend + r4 * np.cos(ang1_aend * np.pi / 180)],
                   y=[by_aend, by_aend - r4 * np.sin(ang1_aend * np.pi / 180)],
                   hoverinfo="skip", line=dict(width=1, dash="dot", color=ext_color), mode="lines"),
        # extension line gamma1 BC X
        go.Scatter(x=[bx_aend, bx_aend - r4 * np.cos(ang2_aend * np.pi / 180)],
                   y=[by_aend, by_aend + r4 * np.sin(ang2_aend * np.pi / 180)],
                   hoverinfo="skip", line=dict(width=1, dash="dot", color=ext_color), mode="lines"),

        # dimensional arc gamma1 X
        go.Scatter(x=[bx_aend + (r4 - radius_offset) * np.sin(beta * np.pi / 180)
                      for beta in range(5, 100 + 1)],
                   y=[by_aend - (r4 - radius_offset) * np.cos(beta * np.pi / 180)
                      for beta in range(5, 100 + 1)],
                   hoverinfo="skip", line=dict(width=1, dash="solid", color=ext_color), mode="lines"),
        # dimensional arc gamma1 arrow left
        go.Scatter(x=[bx_aend + (r4 - radius_offset) * np.cos(ang1_aend * np.pi / 180)],
                   y=[by_aend - (r4 - radius_offset) * np.sin(ang1_aend * np.pi / 180)],
                   hoverinfo="skip", mode="markers",
                   marker=dict(size=7, opacity=1, line=dict(width=0), color=ext_color,
                               symbol="arrow", angle=ang1_aend-180)),
        # dimensional arc gamma1 arrow right
        go.Scatter(x=[bx_aend - (r4 - radius_offset) * np.cos(ang2_aend * np.pi / 180)],
                   y=[by_aend + (r4 - radius_offset) * np.sin(ang2_aend * np.pi / 180)],
                   hoverinfo="skip", mode="markers",
                   marker=dict(size=7, opacity=1, line=dict(width=0), color=ext_color,
                               symbol="arrow", angle=ang2_aend-180)),
        # extension line gamma2 BC X
        go.Scatter(x=[cx_aend, cx_aend + r4 * np.cos(ang2_aend * np.pi / 180)],
                   y=[cy_aend, cy_aend - r4 * np.sin(ang2_aend * np.pi / 180)],
                   hoverinfo="skip", line=dict(width=1, dash="dot", color=ext_color), mode="lines"),
        # extension line gamma2 CD X
        go.Scatter(x=[cx_aend, cx_aend - r4 * np.cos(ang3_aend * np.pi / 180)],
                   y=[cy_aend, cy_aend + r4 * np.sin(ang3_aend * np.pi / 180)],
                   hoverinfo="skip", line=dict(width=1, dash="dot", color=ext_color), mode="lines"),
        # dimensional arc gamma2 X
        go.Scatter(x=[cx_aend - (r4 - radius_offset) * np.cos(beta * np.pi / 180)
                      for beta in range(10, 100 + 1)],
                   y=[cy_aend - (r4 - radius_offset) * np.sin(beta * np.pi / 180)
                      for beta in range(10, 100 + 1)],
                   hoverinfo="skip", line=dict(width=1, dash="solid", color=ext_color), mode="lines"),
        # dimensional arc gamma2 arrow right
        go.Scatter(x=[cx_aend - (r4 - radius_offset) * np.cos(ang3_aend * np.pi / 180)],
                   y=[cy_aend + (r4 - radius_offset) * np.sin(ang3_aend * np.pi / 180)],
                   hoverinfo="skip", mode="markers",
                   marker=dict(size=7, opacity=1, line=dict(width=0), color=ext_color,
                               symbol="arrow", angle=ang3_aend-180)),
        # dimensional arc gamma2 arrow left
        go.Scatter(x=[cx_aend + (r4 - radius_offset) * np.cos(ang2_aend * np.pi / 180)],
                   y=[cy_aend - (r4 - radius_offset) * np.sin(ang2_aend * np.pi / 180)],
                   hoverinfo="skip", mode="markers",
                   marker=dict(size=7, opacity=1, line=dict(width=0), color=ext_color,
                               symbol="arrow", angle=ang2_aend-180)),
        # extension line delta D prim vertical X
        go.Scatter(x=[dx_prim, dx_prim], y=[dy_prim, dy_prim - 20], mode="lines", hoverinfo="skip",
                   line=dict(width=1, dash="dot", color=ext_color)),
        # extension line delta D vertical X
        go.Scatter(x=[dx, dx], y=[dy, dy_prim - 20], mode="lines", hoverinfo="skip",
                   line=dict(width=1, dash="dot", color=ext_color)),
        # extension line delta D prim horizontal X
        go.Scatter(x=[dx_prim, dx + 20], y=[dy_prim, dy_prim], mode="lines", hoverinfo="skip",
                   line=dict(width=1, dash="dot", color=ext_color)),
        # extension line delta D horizontal X
        go.Scatter(x=[dx, dx + 20], y=[dy, dy], mode="lines", hoverinfo="skip",
                   line=dict(width=1, dash="dot", color=ext_color)),
        # dimensional arc alpha bottom X
        go.Scatter(x=curves_alpha[0]["x"], y=curves_alpha[0]["y"], hoverinfo="skip",
                   marker=curves_alpha[0]["marker"],
                   line=dict(width=1, dash="solid", color=ext_color), mode="lines+markers"),
        # dimensional arc alpha top X
        go.Scatter(x=curves_alpha[1]["x"], y=curves_alpha[1]["y"], hoverinfo="skip",
                   marker=curves_alpha[1]["marker"],
                   line=dict(width=1, dash="solid", color=ext_color), mode="lines+markers"),
        # dimensional arc alpha prim bottom X
        go.Scatter(x=curves_alpha_prim[0]["x"], y=curves_alpha_prim[0]["y"], hoverinfo="skip",
                   marker=curves_alpha_prim[0]["marker"],
                   line=dict(width=1, dash="solid", color=ext_color), mode="lines+markers"),
        # dimensional arc alpha prim top X
        go.Scatter(x=curves_alpha_prim[1]["x"], y=curves_alpha_prim[1]["y"], hoverinfo="skip",
                   marker=curves_alpha_prim[1]["marker"],
                   line=dict(width=1, dash="solid", color=ext_color), mode="lines+markers"),
        # dimension line delta D east-north arrow marker
        go.Scatter(x=[dx - arrow_offset], y=[dy - arrow_offset], mode="markers", hoverinfo="skip",
                   marker=dict(size=7, opacity=1, line=dict(width=0), color=ext_color, symbol="arrow",
                               angle=ws_dim_arrow_angle)),
        # dimension line delta D west-south arrow marker
        go.Scatter(x=[dx_prim + arrow_offset], y=[dy_prim + arrow_offset], mode="markers", hoverinfo="skip",
                   marker=dict(size=7, opacity=1, line=dict(width=0), color=ext_color, symbol="arrow",
                               angle=180 + ws_dim_arrow_angle)),
        # dimension line delta D X
        go.Scatter(x=[dx - arrow_offset, dx_prim + arrow_offset], y=[dy - arrow_offset, dy_prim + arrow_offset],
                   mode="lines", line=dict(width=1, dash="solid", color=ext_color), hoverinfo="skip", ),
        # dimension line delta Dx left arrow marker
        go.Scatter(x=[dx], y=[dy_prim - 20], mode="markers", hoverinfo="skip",
                   marker=dict(size=7, opacity=1, line=dict(width=0), color=ext_color, symbol="arrow", angle=90)),
        # dimension line delta Dx right arrow marker
        go.Scatter(x=[dx_prim], y=[dy_prim - 20], mode="markers", hoverinfo="skip",
                   marker=dict(size=7, opacity=1, line=dict(width=0), color=ext_color, symbol="arrow", angle=-90)),
        # dimension line delta Dx X
        go.Scatter(x=[dx, dx_prim], y=[dy_prim - 20, dy_prim - 20], mode="lines", hoverinfo="skip",
                   line=dict(width=1, dash="solid", color=ext_color)),
        # dimension line delta Dy top arrow marker
        go.Scatter(x=[dx + 20], y=[dy], mode="markers", hoverinfo="skip",
                   marker=dict(size=7, opacity=1, line=dict(width=0), color=ext_color, symbol="arrow", angle=0)),
        # dimension line delta Dy bottom arrow marker
        go.Scatter(x=[dx + 20], y=[dy_prim], mode="markers", hoverinfo="skip",
                   marker=dict(size=7, opacity=1, line=dict(width=0), color=ext_color, symbol="arrow", angle=180)),
        # dimension line delta Dy X
        go.Scatter(x=[dx + 20, dx + 20], y=[dy, dy_prim], mode="lines", hoverinfo="skip",
                   line=dict(width=1, dash="solid", color=ext_color)),
        # mechanism alpha zero position
        go.Scatter(x=[0, 40, 75, 100], y=[0, 0, 0, 0], hoverinfo="skip",
                   line=dict(width=4, dash="solid", color="black"), mode="lines"),
        # mechanism alpha position
        go.Scatter(x=[ax, bx, cx, dx], y=[ay, by, cy, dy], hoverinfo="skip",
                   line=dict(width=4, dash="solid", color="black"), mode="lines"),
        # mechanism alpha prim position
        go.Scatter(x=[ax, bx_prim, cx_prim, dx_prim], y=[ay, by_prim, cy_prim, dy_prim], hoverinfo="skip",
                   line=dict(width=4, dash="solid", color="black"), mode="lines"),
        # mechanism alpha almost end position
        go.Scatter(x=[ax, bx_aend, cx_aend, dx_aend], y=[ay, by_aend, cy_aend, dy_aend], hoverinfo="skip",
                   line=dict(width=4, dash="solid", color="black"), mode="lines"),
        # mechanism alpha zero position
        go.Scatter(x=[0, 40, 75, 100], y=[0, 0, 0, 0], hoverinfo="skip",
                   marker=dict(color=mech_color, size=6, line=dict(color="black", width=1)),
                   line=dict(width=2, dash="solid", color=mech_color), mode="markers+lines"),
        # mechanism alpha position
        go.Scatter(x=[ax, bx, cx, dx], y=[ay, by, cy, dy], hoverinfo="skip",
                   marker=dict(color=mech_color, size=6, line=dict(color='black', width=1)),
                   line=dict(width=2, dash="solid", color=mech_color), mode="markers+lines"),
        # mechanism alpha prim position
        go.Scatter(x=[ax, bx_prim, cx_prim, dx_prim], y=[ay, by_prim, cy_prim, dy_prim], hoverinfo="skip",
                   marker=dict(color=mech_color, size=6, line=dict(color='black', width=1)),
                   line=dict(width=2, dash="solid", color=mech_color), mode="markers+lines"),
        # mechanism alpha almost end position
        go.Scatter(x=[ax, bx_aend, cx_aend, dx_aend], y=[ay, by_aend, cy_aend, dy_aend], hoverinfo="skip",
                   marker=dict(color=mech_color, size=6, line=dict(color="black", width=1)),
                   line=dict(width=2, dash="solid", color=mech_color), mode="markers+lines"),
    ]
    layout = go.Layout(
        hovermode=None, hoverlabel=None, showlegend=False, margin=dict(l=5, r=5, t=5, b=5),
        paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(showgrid=False, color="white", showline=False, zeroline=False, visible=False, fixedrange=True),
        yaxis=dict(showgrid=False, color="white", showline=False, zeroline=False, visible=False, fixedrange=True,
                   scaleanchor="x", scaleratio=1, ),
        annotations=[
            # alpha
            dict(showarrow=False, xref="x", yref="y",
                 x=ax - (r2 - radius_offset) * np.cos(ang1 / 2 * np.pi / 180),
                 y=ay + (r2 - radius_offset) * np.sin(ang1 / 2 * np.pi / 180),
                 text=r"$\alpha_n$", textangle=ang1 / 2 - 90, xanchor="center", yanchor="middle",
                 font=dict(color=dim_color, size=11)),
            # alpha prim
            dict(showarrow=False, xref="x", yref="y",
                 x=ax - (r3 - radius_offset) * np.cos(ang2 / 2 * np.pi / 180),
                 y=ay + (r3 - radius_offset) * np.sin(ang2 / 2 * np.pi / 180),
                 text=r"$\alpha_{n+1}$", textangle=ang2 / 2 - 90, xanchor="center", yanchor="middle",
                 font=dict(color=dim_color, size=11)),
            # gamma 1
            dict(showarrow=False, xref="x", yref="y",
                 x=((bx_aend + (r4 - radius_offset) * np.cos(ang1_aend * np.pi / 180)) +
                    (bx_aend - (r4 - radius_offset) * np.cos(ang2_aend * np.pi / 180))) / 2,
                 y=((by_aend - (r4 - radius_offset) * np.sin(ang1_aend * np.pi / 180)) +
                    (by_aend + (r4 - radius_offset) * np.sin(ang2_aend * np.pi / 180))) / 2,
                 text=r"$\gamma_1$", textangle=(ang1_aend + ang2_aend) / 2 - 180, xanchor="center", yanchor="middle",
                 font=dict(color=dim_color, size=11)),
            # gamma 2
            dict(showarrow=False, xref="x", yref="y",
                 x=((cx_aend - (r4 - radius_offset) * np.cos(ang3_aend * np.pi / 180)) +
                    (cx_aend + (r4 - radius_offset) * np.cos(ang2_aend * np.pi / 180))) / 2,
                 y=((cy_aend + (r4 - radius_offset) * np.sin(ang3_aend * np.pi / 180)) +
                    (cy_aend - (r4 - radius_offset) * np.sin(ang2_aend * np.pi / 180))) / 2,
                 text=r"$\gamma_2$", textangle=(ang3_aend + ang2_aend) / 2 - 180, xanchor="center", yanchor="middle",
                 font=dict(color=dim_color, size=11)),
            # D, D prim
            dict(showarrow=False, x=dx, y=dy, xref="x", yref="y", text=r"$ D_n (x_n,y_n)$", xanchor="left",
                 yanchor="bottom", font=dict(color=dim_color, size=11), xshift=5, yshift=5),
            dict(showarrow=False, x=dx_prim, y=dy_prim, xref="x", yref="y", text=r"$ D_{n+1} (x_{n+1},y_{n+1})$",
                 xanchor="right", yanchor="top", font=dict(color=dim_color, size=11), xshift=-5, yshift=-5),
            # A0, B0, C0, D0
            dict(showarrow=False, x=0, y=0, xref="x", yref="y", text=r"$A_0$", xanchor="left", yanchor="bottom",
                 font=dict(color=dim_color, size=11), xshift=5, yshift=5),
            dict(showarrow=False, x=40, y=0, xref="x", yref="y", text=r"$B_0$", xanchor="left", yanchor="bottom",
                 font=dict(color=dim_color, size=11), xshift=5, yshift=5),
            dict(showarrow=False, x=75, y=0, xref="x", yref="y", text=r"$C_0$", xanchor="left", yanchor="bottom",
                 font=dict(color=dim_color, size=11), xshift=5, yshift=5),
            dict(showarrow=False, x=100, y=0, xref="x", yref="y", text=r"$D_0$", xanchor="left", yanchor="bottom",
                 font=dict(color=dim_color, size=11), xshift=5, yshift=5),
            # delta D, delta Dy, delta Dx
            dict(showarrow=False, x=(dx + dx_prim) / 2, y=(dy + dy_prim) / 2, xref="x", yref="y", text=r"$\Delta D$",
                 textangle=-90 + ws_dim_arrow_angle, xanchor="center", yanchor="middle",
                 font=dict(color="#ffd64f", size=8), xshift=-6, yshift=6),
            dict(showarrow=False, x=dx + 20, y=(dy + dy_prim) / 2, xref="x", yref="y", text=r"$\Delta D_y$",
                 textangle=-90, xanchor="right", yanchor="middle", font=dict(color="#ffd64f", size=8),
                 xshift=-2, yshift=0),
            dict(showarrow=False, x=(dx + dx_prim) / 2, y=dy_prim - 20, xref="x", yref="y", text=r"$\Delta D_x$",
                 xanchor="center", yanchor="bottom", font=dict(color="#ffd64f", size=8), xshift=0, yshift=3),
        ]
    )
    figure = go.Figure(data=data, layout=layout)
    return [
        html.Div(id='id24', children=[
            dcc.Graph(id="mechanism-description-chart", config={'displayModeBar': False, 'scrollZoom': False},
                      mathjax=True, figure=figure)]),
    ]


app = dash.Dash(__name__, suppress_callback_exceptions=True)
server = app.server
mathjax = 'https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.4/MathJax.js?config=TeX-MML-AM_CHTML'
app.scripts.append_script({'external_url': mathjax})
app.layout = html.Div([
    dcc.Store(id='memory_container_with_charts'),
    dcc.Store(id='memory_plot2'),
    dcc.Store(id='memory_table_data_and_styles'),
    dcc.Store(id='memory_axis_labels'),
    dcc.Store(id="memory_input_panel_components_values"),
    dcc.Store(id="memory_input_panel_components_styles"),
    dcc.Store(id="memory_model_parameters_trajectory_path"),
    dcc.Store(id="memory_abc"),
    dcc.Store(id="stored_find_parameters_charts"),  # for clientside callback
    dcc.Store(id="stored_mechanism_dimensions_charts"),  # for clientside callback
    dcc.Interval(
        id='interval-component',
        interval=1 * 20,
        n_intervals=0,
        disabled=True
    ),
    html.Div(
        id="id_big-app-container",
        children=[
            html.Div(
                id="id_app-container",
                children=[
                    build_main_tabs(),
                    html.Div(
                        id="id_app-content"
                    )
                ]
            )
        ]
    )
])


def printinfo(func_name, triggered):
    global licznik
    licznik = licznik + 1
    print(f"{licznik}---- EXECUTE {func_name} TRIGGERED BY {triggered}")


@app.callback(
    Output('interval-component', 'n_intervals'),
    Input('interval-component', 'n_intervals')
)
def blablacar(n):
    print(n)


@app.callback(
    Output('id_app-content', 'children'),
    Input('id_app-tabs', 'value'),
    prevent_initial_call=True
)
def render_main_content(tab):
    print("wykonano f-je 'render_main_content'")

    global first_time
    if tab == 'tab0':
        return html.Div([html.H1('Tab content 1')])
    elif tab == 'tab1':
        first_time = True
        init_global_vals()
        return prepare_content_schema_in_tab1()


# @app.callback(
#     Output('show_plot_param1', 'children'),
#     Output('show_plot_param2', 'children'),
#     Output('id_AB', 'disabled'),
#     Output('id_BC', 'disabled'),
#     Output('id_CD', 'disabled'),
#     Output('id_dd1', 'disabled'),
#     Output('id_dd2', 'disabled'),
#     Output('id_dropdown1_div', 'className'),
#     Output('id_dropdown2_div', 'className'),
#     Input('id_btn', 'n_clicks'),
#     Input('show_btn', 'n_clicks'),
#     Input('return_btn', 'n_clicks'),
#     State('id_AB', 'value'),
#     State('id_BC', 'value'),
#     State('id_CD', 'value'),
#     State('id_dd1', 'value'),
#     State('id_dd2', 'value'),
#     State('id_table_show_sols0', 'data'),
#     State('memory_plot1',"data"),
#     State('memory_plot2',"data"),
#     # Input('cont-tabs', 'value'),
#     prevent_initial_call=True
# )
# def calculate_plot(nc1,nc2,nc3, AB, BC, CD, eq_zero_var1, eq_zero_var2, data, oldfig1, oldfig2):
#
#     triggered_id = ctx.triggered_id
#     print('a')
#     print(triggered_id)
#     if triggered_id == 'id_btn':
#         built_plot1, built_plot2 = build_find_parameters_plot(nc1,AB,BC,CD,eq_zero_var1, eq_zero_var2)
#         return built_plot1, built_plot2, no_update, no_update, no_update, no_update, no_update, no_update, no_update
#     elif triggered_id == 'show_btn':
#         if nc2 != 0:
#             for param in table_header_names:
#                 program = param + "=data[0]['" + param + "']"
#                 exec(program, locals(), globals())
#             joint_trajectory = calculate_model(AB, BC, CD, mx, my, nx, ny, px, py, sx, sy)
#             built_plot3, built_plot4 = build_mechanism_plot(nc2,joint_trajectory, AB, BC, CD)
#             return built_plot3, built_plot4, True, True, True, True, True, "dropdown-div-disabled", "dropdown-div-disabled"
#     elif triggered_id == 'return_btn':
#         return oldfig1, oldfig2, False, False, False, False, False, "dropdown-div-enabled", "dropdown-div-enabled"
#
#
#
# @app.callback(
#     Output('id_table_show_sols0', 'data'), # ok
#     Output("id_graph1", "clickData"), # not ok
#     Output("id_graph2", "clickData"), # not ok
#     Output("show_btn", "disabled"), # ok
#     Output("id_graph1", "figure"), # not ok
#     Output("id_graph2", "figure"), # not ok
#     Output("return_btn", "disabled"), # ok
#     Output("id_btn", "disabled"), # ok
#     Input('id_btn', 'n_clicks'), # ok
#     State("id_graph1", "figure"), # not ok
#     State("id_graph2", "figure"), # not ok
#     State('id_table_show_sols0', 'data'), # ok
#     Input("id_graph1", "clickData"), # not ok
#     Input("id_graph2", "clickData"), # not ok
#     Input('id_dd1', 'value'), # ok
#     Input('id_dd2', 'value'), # ok
#     Input('id_AB', 'value'), # ok
#     Input('id_BC', 'value'), # ok
#     Input('id_CD', 'value'), # ok
#     Input('show_btn', 'n_clicks'), # ok
#     Input('return_btn', 'n_clicks'), # ok
#     Input('id_btn', 'disabled'), # ok
#     Input('cont-tabs', 'value')
#     # prevent_initial_call=True
# )
# def show_table(nc3, fig1, fig2, old_data, clickData1, clickData2, d1val, d2val, AB, BC, CD, nc1, nc2, calc_btn_dis, tabname):
#     triggered_id = ctx.triggered_id
#     # triggered_prop_id = ctx.triggered_prop_ids
#     global show_button_release3, show_button_release1, show_button_release2, first_time
#
#     empty_data = [dict(**{param: '-' for param in table_header_names})]
#     print('b')
#     print(triggered_id)
#     if triggered_id == "cont-tabs":
#         if tabname == 'dimensional_drawing':
#             return no_update, no_update, no_update, True, no_update, no_update, True, True
#         elif tabname == 'choose_solutions':
#             return no_update, no_update, no_update, True, no_update, no_update, True, False
#         elif tabname == 'model_parameters':
#             return no_update, no_update, no_update, True, no_update, no_update, True, True
#     elif triggered_id == 'show_btn':
#         if nc1 != 0:
#             return no_update, no_update, no_update, True, no_update, no_update, False, True
#     elif triggered_id == "id_btn":
#         if nc3 != 0:
#             return no_update, no_update, no_update, no_update, no_update, no_update, no_update, True
#     elif (triggered_id == 'id_AB') or (triggered_id == 'id_BC') or (triggered_id == 'id_CD'):
#         return empty_data, None, None, True, no_update, no_update, no_update, False
#     elif (triggered_id == 'id_dd1') or (triggered_id == 'id_dd2'):
#         show_button_release3 = 0
#         show_button_release1 = 0
#         show_button_release2 = 0
#         disabled_btn = True
#         return empty_data, no_update, no_update, disabled_btn, no_update, no_update, no_update, False
#
#     # Jeżeli wykreślono scatter3d
#     elif show_button_release3 == 1:
#         cd1 = clickData1
#         cd2 = clickData2
#         x1_label, y1_label, z1_label, zero1_label, x2_label, y2_label, z2_label, zero2_label = get_labels(d1val, d2val)
#
#         if ((clickData1 is None) and (clickData2 is None)) or (calc_btn_dis == False):
#             new_data = empty_data
#         else:
#             # uzupełnienie 1 czesci tabeli
#             if clickData1 is not None:
#                 new_data = update_table(old_data, clickData1, x1_label, y1_label, z1_label, d1val)
#                 show_button_release1 = 1
#             # uzupełnienie 2 częsci tabeli
#             if clickData2 is not None:
#                 new_data = update_table(old_data, clickData2, x2_label, y2_label, z2_label, d2val)
#                 show_button_release2 = 1
#
#         # Jeżeli tabelka jest wypełniona to odblokuj przycisk SHOW    ?????????????????????????????????????
#         if (show_button_release2 != 0) and (show_button_release1 != 0) and (clickData1 is not None) and (clickData2 is not None) and (calc_btn_dis == True):
#             disabled_btn = False
#         else:
#             disabled_btn = True
#
#         if first_time == True:
#             first_time = False
#             return no_update, no_update, no_update, no_update, no_update, no_update, no_update, True
#
#         elif (triggered_id == "id_graph1") and (calc_btn_dis == True):
#             if cd1 is not None:
#                 old_clist1 = fig1['data'][0]['marker']['color']
#                 old_clist2 = fig2['data'][0]['marker']['color']
#                 color_list1, color_list2 = color_plot(1,cd1,d1val, d2val, AB, BC, CD, old_clist1, old_clist2)
#                 fig1['data'][0]['marker']['color'] = color_list1
#                 fig2['data'][0]['marker']['color'] = color_list2
#
#         elif (triggered_id == "id_graph2") and (calc_btn_dis == True):
#             if cd2 is not None:
#                 old_clist1 = fig1['data'][0]['marker']['color']
#                 old_clist2 = fig2['data'][0]['marker']['color']
#                 color_list1, color_list2 = color_plot(2, cd2, d1val, d2val, AB, BC, CD, old_clist1, old_clist2)
#                 fig2['data'][0]['marker']['color'] = color_list2
#                 fig1['data'][0]['marker']['color'] = color_list1
#         return new_data, no_update, no_update, disabled_btn, fig1, fig2, True, no_update

#
# @app.callback(
#     Output('memory_plot1','data'),
#     Output('memory_plot2','data'),
#     State('show_plot_param1', 'children'),
#     State('show_plot_param2', 'children'),
#     Input("id_graph1", "clickData"),
#     Input("id_graph2", "clickData"),
# )
# def storage_figures(fig1,fig2,cd1,cd2):
#     triggered_id = ctx.triggered_id
#
#     if (triggered_id == 'id_graph1') and (cd1 is not None):
#         # print(f"ZAPISANO FIGURY W STORAGE WYWOŁUJĄC  {triggered_id}")
#         return fig1, fig2
#     elif (triggered_id == 'id_graph2') and (cd2 is not None):
#         # print(f"ZAPISANO FIGURY W STORAGE WYWOŁUJĄC  {triggered_id}")
#         return fig1, fig2
#     else:
#         return no_update, no_update


# @app.callback(
#     Output('cont2',"children"),
#     Input('id_btn','n_clicks'),
#     Input('id_AB', 'value')
# )
# def lolizka(nc,val):
#     triggered_id = ctx.triggered_id
#     print('c')
#     print(triggered_id)
#     return no_update


# @app.callback(
#     Output("id_app-container", 'children'),
#     # Output('id_btn', 'disabled'),
#     Input('id_graph3', 'figure'),
#     # State('id_graph3', 'layout'),
#     # State('return_btn','disabled'),
#     # prevent_initial_call=True
# )
# def lols(x):
#     # hD["points"][0]["marker.color"] = 'red'
#     # hD["points"][0]["marker.size"] = 8
#     print(x)
#     return no_update
# @app.callback(
#     Output("memory_model_parameters_trajectory_path", "data"),
#     Input('cont-tab3', 'disabled'),
#     State('id_AB', 'value'),
#     State('id_BC', 'value'),
#     State('id_CD', 'value'),
#     State('memory_table_data_and_styles', 'data'),
# )
# def store_model_parameters_trajectory_path(disabled_tab, ab, bc, cd, mtdas):
#     if disabled_tab is False:
#         printinfo("store_model_parameters_trajectory_path", ctx.triggered_id)
#         return calculate_model(ab, bc, cd, mtdas["data"][0]).to_json()
#     else:
#         return no_update
@app.callback(
    Output("slider", "value"),
    Output('interval-component', 'disabled'),
    State("slider", "value"),
    Input('interval-component', 'n_intervals'),
    Input("button_model_parameters_start", "n_clicks"),
    Input("button_model_parameters_pause", "n_clicks"),
    Input("button_model_parameters_next", "n_clicks"),
    prevent_initial_call=True)
def update_slider_value_model_parameters(old_slider_value, n_intervals, start_nclicks, pause_nclicks, next_nclicks):
    print(ctx.triggered_id)
    step = 0.1
    if old_slider_value == 90:
        print("A")
        return no_update, True
    elif ctx.triggered_id == "button_model_parameters_start" and start_nclicks:
        print("B")
        return old_slider_value + step, False
    elif ctx.triggered_id == "button_model_parameters_pause" and pause_nclicks:
        print("C")
        return no_update, True
    elif ctx.triggered_id == "button_model_parameters_next" and next_nclicks:
        print("D")
        return old_slider_value + step, True
    elif ctx.triggered_id == "interval-component":
        print("E")
        return old_slider_value + step, no_update
    else:
        return no_update, no_update


@app.callback(
    Output("mechanism-positions-plot", "figure"),
    Output("mechanism-vectors-plot", "figure"),
    Output("id_graph5", "figure"),
    Output("id_graph6", "figure"),
    Output("id_graph11", "figure"),
    Output("mechanism-angles-dynamic", "figure"),
    Input("slider", "value"),
    State("memory_model_parameters_trajectory_path", "data"),
    State("mechanism-positions-plot", "figure"),
    State("mechanism-vectors-plot", "figure"),
    State("id_graph5", "figure"),
    State("id_graph6", "figure"),
    State("id_graph11", "figure"),
    State("mechanism-angles-dynamic", "figure")
)
def update_model_parameters_plot(slider_value, mmptp, figure1, figure2, figure3, figure4, figure5, figure6):
    # The function updates charts in the "model parameters" tab based on calls via slider
    # (check: "update_slider_value_model_parameters")
    # "si" is the number of current sample. Slider value step is 0.1 degree, so "si" is multiplied by 10
    step = 0.1
    multiplexer = 1 / step
    si = int(slider_value * multiplexer)

    def update_mechanism_angles_dynamic_chart():
        figure6['data'][0]['y'] = [90 + mmptp['kABC'][si], 90 + mmptp['kBCD'][si]]
        figure6['layout']['title']['text'] = r"$\alpha\text{ = }" + f"{mmptp['alpha'][si]:.1f}" + r"^\circ$"
        figure6['data'][0]['text'] = [f"{mmptp['kABC'][si]:.2f}\u00B0", f"{mmptp['kBCD'][si]:.2f}\u00B0"]
        return figure6

    def update_bar_charts():
        # The function updates three bar charts: "x-direction velocity", "y-direction velocity" and "resultant velocity"
        # "dBx" is difference between values of x-position in steps "n" and "n-1"
        # "mv" prefix - "modified value" <- prepared value to show bars with appropriate range
        #   (calculated by dBx/max(dBx), and show as value <0,1>)
        # "tv" prefix - "text value" <- prepared text label with original value
        key_names_list = [["dBx", "dCx", "dDx", "dFx", "dGx", "dHx"],
                          ["dBy", "dCy", "dDy", "dFy", "dGy", "dHy"],
                          ["dB", "dC", "dD", "dF", "dG", "dH"]]
        bar_charts_list = []
        for idx, fig in enumerate([figure3, figure4, figure5]):
            key_names = key_names_list[idx]
            fig['data'][0]['y'] = [mmptp["mv" + key_names[0]][si], mmptp["mv" + key_names[1]][si],
                                   mmptp["mv" + key_names[2]][si], mmptp["mv" + key_names[3]][si],
                                   mmptp["mv" + key_names[4]][si], mmptp["mv" + key_names[5]][si]]
            fig['data'][0]['text'] = [mmptp["tv" + key_names[0]][si], mmptp["tv" + key_names[1]][si],
                                      mmptp["tv" + key_names[2]][si], mmptp["tv" + key_names[3]][si],
                                      mmptp["tv" + key_names[4]][si], mmptp["tv" + key_names[5]][si]]
            bar_charts_list.append(fig)
        return bar_charts_list

    def update_position(fig):
        # Auxiliary function; function is called when "position chart" or "vector chart" is updated.
        # The function updates traces: ABG, FBC, HCD, EF, GH
        key_list_names = [[["Ax", "Bx", "Gx"], ["Ay", "By", "Gy"]], [["Fx", "Bx", "Cx"], ["Fy", "By", "Cy"]],
                          [["Hx", "Cx", "Dx"], ["Hy", "Cy", "Dy"]], [["Gx", "Hx"], ["Gy", "Hy"]],
                          [["Ex", "Fx"], ["Ey", "Fy"]]]
        for trace_no in range(0, 5):
            for idx, ax_name in enumerate(["x", "y"]):
                fig["data"][trace_no][ax_name] = [mmptp[element][si] for element in key_list_names[trace_no][idx]]
        return fig

    def update_position_chart():
        # The function updates "positional chart". Only segments ABG, FBC, HCD, EF, GH are updated.
        return update_position(figure1)

    def update_vector_chart():
        # The function updates "vector chart". Updated are segments and vectors.
        fig = update_position(figure2)
        # x_list = [[mmptp["Bx"][si], mmptp["vhBx"][si]],
        #           [mmptp["Cx"][si], mmptp["vhCx"][si]],
        #           [mmptp["Dx"][si], mmptp["vhDx"][si]],
        #           [mmptp["Fx"][si], mmptp["vhFx"][si]],
        #           [mmptp["Gx"][si], mmptp["vhGx"][si]],
        #           [mmptp["Hx"][si], mmptp["vhHx"][si]],
        #           [mmptp["Bx"][si], mmptp["vvBx"][si]],
        #           [mmptp["Cx"][si], mmptp["vvCx"][si]],
        #           [mmptp["Dx"][si], mmptp["vvDx"][si]],
        #           [mmptp["Fx"][si], mmptp["vvFx"][si]],
        #           [mmptp["Gx"][si], mmptp["vvGx"][si]],
        #           [mmptp["Hx"][si], mmptp["vvHx"][si]]]
        #
        # y_list = [[mmptp["By"][si], mmptp["vhBy"][si]],
        #           [mmptp["Cy"][si], mmptp["vhCy"][si]],
        #           [mmptp["Dy"][si], mmptp["vhDy"][si]],
        #           [mmptp["Fy"][si], mmptp["vhFy"][si]],
        #           [mmptp["Gy"][si], mmptp["vhGy"][si]],
        #           [mmptp["Hy"][si], mmptp["vhHy"][si]],
        #           [mmptp["By"][si], mmptp["vvBy"][si]],
        #           [mmptp["Cy"][si], mmptp["vvCy"][si]],
        #           [mmptp["Dy"][si], mmptp["vvDy"][si]],
        #           [mmptp["Fy"][si], mmptp["vvFy"][si]],
        #           [mmptp["Gy"][si], mmptp["vvGy"][si]],
        #           [mmptp["Hy"][si], mmptp["vvHy"][si]]]

        x_list = [[mmptp["Bx"][si], mmptp["vhBx"][si]],
                  [mmptp["Bx"][si], mmptp["vvBx"][si]],
                  [mmptp["Cx"][si], mmptp["vhCx"][si]],
                  [mmptp["Cx"][si], mmptp["vvCx"][si]],
                  [mmptp["Dx"][si], mmptp["vhDx"][si]],
                  [mmptp["Dx"][si], mmptp["vvDx"][si]],
                  [mmptp["Fx"][si], mmptp["vhFx"][si]],
                  [mmptp["Fx"][si], mmptp["vvFx"][si]],
                  [mmptp["Gx"][si], mmptp["vhGx"][si]],
                  [mmptp["Gx"][si], mmptp["vvGx"][si]],
                  [mmptp["Hx"][si], mmptp["vhHx"][si]],
                  [mmptp["Hx"][si], mmptp["vvHx"][si]]]

        y_list = [[mmptp["By"][si], mmptp["vhBy"][si]],
                  [mmptp["By"][si], mmptp["vvBy"][si]],
                  [mmptp["Cy"][si], mmptp["vhCy"][si]],
                  [mmptp["Cy"][si], mmptp["vvCy"][si]],
                  [mmptp["Dy"][si], mmptp["vhDy"][si]],
                  [mmptp["Dy"][si], mmptp["vvDy"][si]],
                  [mmptp["Fy"][si], mmptp["vhFy"][si]],
                  [mmptp["Fy"][si], mmptp["vvFy"][si]],
                  [mmptp["Gy"][si], mmptp["vhGy"][si]],
                  [mmptp["Gy"][si], mmptp["vvGy"][si]],
                  [mmptp["Hy"][si], mmptp["vhHy"][si]],
                  [mmptp["Hy"][si], mmptp["vvHy"][si]]]

        marker_symbol_list = [
            mmptp["vhsB"][si], mmptp["vvsB"][si],
            mmptp["vhsC"][si], mmptp["vvsC"][si],
            mmptp["vhsD"][si], mmptp["vvsD"][si],
            mmptp["vhsF"][si], mmptp["vvsF"][si],
            mmptp["vhsG"][si], mmptp["vvsG"][si],
            mmptp["vhsH"][si], mmptp["vvsH"][si],
        ]

        for idx, trace_no in enumerate(range(-12, 0, 1)):
            fig['data'][trace_no]['x'] = x_list[idx]
            fig['data'][trace_no]['y'] = y_list[idx]
            fig['data'][trace_no]['marker']['symbol'] = marker_symbol_list[trace_no]

        return fig

    return [update_position_chart(), update_vector_chart()] + update_bar_charts() + \
        [update_mechanism_angles_dynamic_chart()]


@app.callback(
    Output('cont3', 'children'), Output('id_AB', 'disabled'), Output('id_BC', 'disabled'), Output('id_CD', 'disabled'),
    Output('id_dd1', 'disabled'), Output('id_dd2', 'disabled'),
    Output('id_dropdown1_div', 'className'), Output('id_dropdown2_div', 'className'),
    Output('id_btn', 'disabled'), Output('cont-tab3', 'disabled'),
    Output("id_table_show_sols0", 'style_cell'),
    Output("id_table_show_sols0", 'style_header'),
    Output("id_table_show_sols0", 'data'),
    Output("memory_model_parameters_trajectory_path", "data"),
    Input('cont-tabs', 'value'), Input('id_btn', 'n_clicks'), State('id_btn', 'disabled'),
    Input('id_dd1', 'value'), Input('id_dd2', 'value'), Input('id_AB', 'value'),
    Input('id_BC', 'value'), Input('id_CD', 'value'),
    State('memory_container_with_charts', 'data'),
    Input('memory_table_data_and_styles', 'data'),
    State("memory_input_panel_components_styles", "data"),
    prevent_initial_call=True
)
def manage_content_style_in_tab1(sub_tab_name, _nc_calc, _btn_calc_disabled, _dd1, _dd2, ab, bc, cd, mcwc, mtdas,
                                 mipcs):
    printinfo("manage_content_style_in_tab1", ctx.triggered_id)

    def update_content_after_change_input_parameters():
        return no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, \
            False, no_update, tabstyle_cell_disabled, tabstyle_header_disabled, empty_table_data, no_update

    def update_content_after_click_calc_button():
        return no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, \
            True, no_update, tabstyle_cell_enabled, tabstyle_header_enabled, empty_table_data, no_update

    def show_content_after_change_tab():
        return prepare_dimensional_drawing_content(), no_update, no_update, no_update, no_update, no_update, \
            no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update

    def update_content_after_change_data_table():
        if "-" in mtdas["data"][0].values():
            return no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, \
                no_update, no_update, no_update, mtdas['data'], no_update
        else:
            return no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, \
                False, no_update, no_update, mtdas['data'], no_update

    def show_dimensional_drawing_content():
        return prepare_dimensional_drawing_content(), True, True, True, True, True, "dropdown-div-disabled", \
            "dropdown-div-disabled", True, no_update, tabstyle_cell_disabled, tabstyle_header_disabled, no_update, \
            no_update

    def show_saved_choose_solutions_content():
        print(mcwc)
        return [mcwc] + list(mipcs.values()) + [no_update]

    def show_basic_choose_solutions_content():
        return prepare_choose_solutions_content(), False, False, False, False, False, "dropdown-div-enabled", \
            "dropdown-div-enabled", False, no_update, tabstyle_cell_disabled, tabstyle_header_disabled, \
            empty_table_data, no_update

    def show_model_parameters_content():
        jt_as_dict, div_structure = prepare_model_parameters_content()
        return div_structure, True, True, True, True, True, "dropdown-div-disabled", \
            "dropdown-div-disabled", True, no_update, tabstyle_cell_disabled, tabstyle_header_disabled, no_update, \
            jt_as_dict

    def prepare_dimensional_drawing_content():
        return [html.Div(id="id_pct_dim1", children=[html.Img(id="id_pct_dim2", src=app.get_asset_url("wymiary.png"))])]

    def prepare_choose_solutions_content():
        return [
            html.Div(id="show_sol_cont", children=[html.Div(id='show_plot_param1'), html.Div(id='show_plot_param2')]),
            html.Div(id="show_sol_cont3",
                     children=[html.Div(id='show_plot_param13'), html.Div(id='show_plot_param23')])]

    def prepare_model_parameters_content():
        joint_trajectory = calculate_model(ab, bc, cd, mtdas["data"][0])

        div_structure = [
            html.Div(id="show_sol_cont44", children=[
                html.Div(id="slider_container", children=build_slider_mechanism_params_content(joint_trajectory)),
                html.Div(id="id3", children=[build_mechanism_position_chart("normal", joint_trajectory, ab, bc, cd)]),
                html.Div(id="id4", children=[build_mechanism_position_chart("vectors", joint_trajectory, ab, bc, cd)])
            ]),
            html.Div(id="show_sol_cont77", children=[
                html.Div(id="container_header", children=build_description_chart()),
                html.Div(id="container_1_line", children=[
                    html.Div(id="container_1_line_1_element", children=build_static_velocity_chart(joint_trajectory)),
                    html.Div(id="container_1_line_2_element", children=build_static_angles_chart(joint_trajectory)),

                ]),
                html.Div(id="container_2_line", children=[
                    html.Div(id="container_2_line_1_element", children=[
                        build_bar_chart("x", joint_trajectory)
                    ]),
                    html.Div(id="container_2_line_2_element", children=[
                        build_bar_chart("y", joint_trajectory)
                    ]),
                    html.Div(id="container_2_line_3_element", children=[
                        build_bar_chart("resultant", joint_trajectory)
                    ])
                ])

            ])
        ]
        return joint_trajectory, div_structure

    triggered_id = ctx.triggered_id
    empty_table_data = [dict(**{name: '-' for name in table_header_names})]
    if triggered_id is None:
        return show_content_after_change_tab()
    elif triggered_id == 'cont-tabs':
        if sub_tab_name == 'dimensional_drawing':
            return show_dimensional_drawing_content()
        elif sub_tab_name == 'choose_solutions':
            if mcwc is not None:
                print("code a")
                return show_saved_choose_solutions_content()
            else:
                return show_basic_choose_solutions_content()
        elif sub_tab_name == 'model_parameters':
            return show_model_parameters_content()
    elif triggered_id == 'memory_table_data_and_styles' and mtdas is not None:
        return update_content_after_change_data_table()
    elif triggered_id == 'id_btn':
        return update_content_after_click_calc_button()
    elif triggered_id in ["id_AB", "id_BC", "id_CD", "id_dd1", "id_dd2"]:
        return update_content_after_change_input_parameters()


@app.callback(
    # Output('show_plot_param1', 'children'), Output('show_plot_param2', 'children'),
    # Output('show_plot_param13', 'children'), Output('show_plot_param23', 'children'),
    Output("stored_find_parameters_charts", "data"),
    Output("stored_mechanism_dimensions_charts", "data"),
    Input('id_btn', 'n_clicks'), State('id_dd1', 'value'), State('id_dd2', 'value'),
    State('id_AB', 'value'), State('id_BC', 'value'), State('id_CD', 'value'),
    prevent_initial_call=True
)
def build_charts_on_choose_solutions_content(nc_calc, dd1, dd2, ab, bc, cd):
    printinfo("build_charts_on_choose_solutions_content", ctx.triggered_id)

    if nc_calc != 0:
        return build_find_parameters_plot(ab, bc, cd, dd1, dd2), \
            build_mechanism_dimensions_chart(ab, bc, cd, dd1, dd2, None)

# @app.callback(
#     Output('show_plot_param1', 'children'),
#     Output('show_plot_param2', 'children'),
#     Input("stored_find_parameters_charts", "data"),
#     prevent_initial_call=True
#
# )
# def build_choose_solutions_charts_by_stored_data(figures):
#     print("build_choose_solutions_charts_by_stored_data")
#
#     return (dcc.Graph(clear_on_unhover=True, config={'displayModeBar': False, 'scrollZoom': False}, id='id_graph1',
#                       figure=figures[0]),
#             dcc.Graph(clear_on_unhover=True, config={'displayModeBar': False, 'scrollZoom': False}, id='id_graph2',
#                       figure=figures[1]))

# @app.callback(
#     Output('show_plot_param13', 'children'),
#     Output('show_plot_param23', 'children'),
#     Input("stored_mechanism_dimensions_charts", "data"),
#     prevent_initial_call=True
#
# )
# def build_mechanism_dimensions_charts_by_stored_data(figures):
#     print("build_mechanism_dimensions_charts_by_stored_data")
#
#     return figures[0], figures[1]


# CLIENTSIDE CALLBACK FOR SHOW 'MECHANISM DIMENSIONS CHARTS'
app.clientside_callback(
    """
    function(data) {
        return [data[0], data[1]];
    }
    """,
    [Output('show_plot_param13', 'children'),
     Output('show_plot_param23', 'children')],
    Input("stored_mechanism_dimensions_charts", "data"),
)

# CLIENTSIDE CALLBACK FOR SHOW 'FIND PARAMETERS CHARTS'
app.clientside_callback(
    """
    function(data) {
        return [data[0], data[1]];
    }
    """,
    [Output('show_plot_param1', 'children'),
     Output('show_plot_param2', 'children')],
    Input("stored_find_parameters_charts", "data"),
)

# HOVER POINTS
app.clientside_callback(
    """
    function(hd1, hd2, fig1, fig2) {
    
        data_length = fig1.data[0].marker.color.length;
        const new_data = fig1.data;
        const default_color_array = Array(data_length).fill("#85c5ed");
        const default_size_array = Array(data_length).fill(4);
        
        if (typeof hd1 === "undefined") {
            // set default color array
            new_data[0].marker.color = default_color_array;
            
            // set default size array
            new_data[0].marker.size = default_size_array;
            
            // update figure
            new_fig = Object.assign({}, fig1, {'data': new_data});
            return [new_fig, fig2];
        } else {
            // prepare color array
            new_color_array = default_color_array;
            new_color_array[hd1.points[0].pointNumber] = "#1dcfb1";
            
            // prepare size array
            new_size_array = default_size_array;
            new_size_array[hd1.points[0].pointNumber] = 10;
            
            // set new color and size arrays
            new_data[0].marker.color = new_color_array;
            new_data[0].marker.size = new_size_array;
            
            // update figure
            new_fig = Object.assign({}, fig1, {'data': new_data});
            return [new_fig, fig2];
            
        }

    }
    """,
    [Output('id_graph1', 'figure'),
     Output('id_graph2', 'figure')],
    Input("id_graph1", "hoverData"),
    Input("id_graph2", "hoverData"),
    State("id_graph1", "figure"),
    State("id_graph2", "figure"),
    # prevent_initial_call=True

)


# @app.callback(
#     Output("id_graph1", "figure"), Output("id_graph2", "figure"),
#     Output("id_graph1", "clickData"), Output("id_graph2", "clickData"),
#     Input("id_graph1", "clickData"), Input("id_graph2", "clickData"),
#     State("id_graph1", "figure"), State("id_graph2", "figure"),
#     Input('memory_axis_labels', 'data'), Input("id_btn", "disabled"),
#     Input("id_graph1", "hoverData"), Input("id_graph2", "hoverData"),
#     prevent_initial_call=True
# )
# def select_active_points(cd1, cd2, fig1, fig2, labels, calc_btn_disabled, hd1, hd2):
#     printinfo("select_active_points", ctx.triggered_id)
#
#     triggered_id = ctx.triggered_id
#
#     # -- clear graph when user changes the input parameters
#     if calc_btn_disabled == False:
#         fig1['data'][0]['x'], fig1['data'][0]['y'], fig1['data'][0]['z'] = [], [], []
#         fig2['data'][0]['x'], fig2['data'][0]['y'], fig2['data'][0]['z'] = [], [], []
#         return fig1, fig2, None, None
#
#     # -- change markers style to default when inputs are changed
#     if triggered_id == "memory_input_panel_components_values":
#         fig1['data'][0]['marker']['color'] = \
#             [marker_styles["FindParamsScatter3d"]["default_marker_color"] for i in fig1['data'][0]['marker']['color']]
#         fig1['data'][0]['marker']['size'] = \
#             [marker_styles["FindParamsScatter3d"]["default_marker_size"] for i in fig1['data'][0]['marker']['size']]
#         fig2['data'][0]['marker']['color'] = \
#             [marker_styles["FindParamsScatter3d"]["default_marker_color"] for i in fig2['data'][0]['marker']['color']]
#         fig2['data'][0]['marker']['size'] = \
#             [marker_styles["FindParamsScatter3d"]["default_marker_size"] for i in fig2['data'][0]['marker']['size']]
#         return fig1, fig2, None, None
#
#     # -- change markers style when points are clicked or hovered
#     elif list(ctx.triggered_prop_ids.keys())[0] in ['id_graph1.hoverData', 'id_graph2.hoverData',
#                                                     'id_graph1.clickData', 'id_graph2.clickData']:
#         # get coordinates of points necessary to create an instance of the "Marker" class
#         x1, y1, z1 = fig1["data"][0]["x"], fig1["data"][0]["y"], fig1["data"][0]["z"]
#         x2, y2, z2 = fig2["data"][0]["x"], fig2["data"][0]["y"], fig2["data"][0]["z"]
#         # prepare labels lists
#         labels_list1 = [labels["x1l"], labels["y1l"], labels["z1l"], labels["zero1l"]]
#         labels_list2 = [labels["x2l"], labels["y2l"], labels["z2l"], labels["zero2l"]]
#
#         # create an instance of the "Marker" class and call "update" to update markers style
#         # id 'update' method first cd and labels_list are own, second are foreign
#         markers_fig1 = Markers(fig1["data"][0]["marker"], x1, y1, z1)
#         m1 = markers_fig1.update(hd1, cd1, cd2, labels_list1, labels_list2, marker_styles["FindParamsScatter3d"])
#         markers_fig2 = Markers(fig2["data"][0]["marker"], x2, y2, z2)
#         m2 = markers_fig2.update(hd2, cd2, cd1, labels_list2, labels_list1, marker_styles["FindParamsScatter3d"])
#         fig1["data"][0]["marker"] = m1
#         fig2["data"][0]["marker"] = m2
#         return fig1, fig2, no_update, no_update
#     else:
#         return no_update, no_update, no_update, no_update


# @app.callback(
#     Output("id_graph7", "figure"), Output("id_graph8", "figure"),
#     Input("id_table_show_sols0", 'data'), State('id_dd1', 'value'), State('id_dd2', 'value'), State('id_AB', 'value'),
#     State('id_BC', 'value'), State('id_CD', 'value'), Input("id_graph1", "hoverData"), Input("id_graph2", "hoverData"),
#     State('memory_axis_labels', 'data'), Input("id_graph7", "hoverData"), Input("id_graph8", "hoverData"),
#     State("id_graph7", "figure"), State("id_graph8", "figure"),
#     prevent_initial_call=True)
# def update_mech_dims_chart_by_hovered_points(
#         table_data, equal_zero_variable1_name_input, equal_zero_variable2_name_input, ab_input, bc_input, cd_input, hd1,
#         hd2, labels, hd7, hd8, mechanism_abc_chart, mechanism_bcd_chart):
#     printinfo("update_mech_dims_chart_by_hovered_points", ctx.triggered_id)
#
#     # Change markers style in 'id_graph7' or 'id_graph8' when point on 'id_graph7' or 'id_graph8' is hovered.
#     # "Try" avoids an error when the chart is loaded for the first time
#     if list(ctx.triggered_prop_ids.keys())[0] in ['id_graph7.hoverData', 'id_graph8.hoverData']:
#         try:
#             markers_fig_abc = Markers(mechanism_abc_chart["data"][-1]["marker"])
#             m1 = markers_fig_abc.update(hd=hd7, marker_styles=marker_styles["MechanismDimensions"])
#             mechanism_abc_chart["data"][-1]["marker"] = m1
#             markers_fig_bcd = Markers(mechanism_bcd_chart["data"][-1]["marker"])
#             m2 = markers_fig_bcd.update(hd=hd8, marker_styles=marker_styles["MechanismDimensions"])
#             mechanism_bcd_chart["data"][-1]["marker"] = m2
#             return mechanism_abc_chart, mechanism_bcd_chart
#         except:
#             return no_update, no_update
#
#     # If function is called by 'id_graph1' or 'id_graph2' and any of points is hovered then update the table data
#     # and return 'build_mechanism_dimensions_chart()' with modified table data.
#     # Otherwise, return 'build_mechanism_dimensions_chart()' with not modified table data
#     elif list(ctx.triggered_prop_ids.keys())[0] in ['id_graph1.hoverData', 'id_graph2.hoverData'] and any([hd1, hd2]):
#         if hd1 is not None:
#             xl, yl, zl, zerol = labels["x1l"], labels["y1l"], labels["z1l"], labels["zero1l"]
#             hd = hd1
#         else:
#             xl, yl, zl, zerol = labels["x2l"], labels["y2l"], labels["z2l"], labels["zero2l"]
#             hd = hd2
#
#         table_data[0][xl] = hd['points'][0]['x']
#         table_data[0][yl] = hd['points'][0]['y']
#         table_data[0][zl] = hd['points'][0]['z']
#         table_data[0][zerol] = 0
#
#     return build_mechanism_dimensions_chart(
#         ab_input, bc_input, cd_input, equal_zero_variable1_name_input, equal_zero_variable2_name_input, table_data)


@app.callback(
    Output('memory_container_with_charts', 'data'),
    Input('id_graph1', 'figure'), Input('id_graph2', 'figure'), State('cont3', 'children'),
    State('cont-tab3', 'disabled'))
def store_container_with_charts(_fig1, _fig2, child, _model_parameters_tab_disabled):
    printinfo("store_container_with_charts", ctx.triggered_id)
    # Store container with charts 'id_graph1', 'id_graph2', 'id_graph7', 'id_graph8'
    # Stored data is used to save a part of view. The view will load when tab "tab1" are switched again
    # (check: 'manage_content_style_in_tab1')
    return child


@app.callback(
    Output('memory_axis_labels', 'data'),
    Input('id_btn', 'n_clicks'), State('id_dd1', 'value'), State('id_dd2', 'value'),
    prevent_initial_call=True)
def store_axis_labels(_, eq_zero_var1, eq_zero_var2):
    printinfo("store_axis_labels", ctx.triggered_id)

    # Stored axis labels from solutions charts
    # Stored data is used to:
    #  - markers style change in 'id_graph1' and 'id_graph2' after click the point on the graph
    #  (check: 'select_active_points')
    #  - mechanism charts update when hovered solutions charts (check: 'update_mech_dims_chart_by_hovered_points')
    def get_labels(eq_zero_var):
        x_label, y_label, z_label, zero_label = None, None, None, None
        if eq_zero_var in ["mx", "my", "nx", "ny"]:
            if eq_zero_var == "mx":
                x_label, y_label, z_label, zero_label = "nx", "ny", "my", "mx"
            elif eq_zero_var == "my":
                x_label, y_label, z_label, zero_label = "nx", "ny", "mx", "my"
            elif eq_zero_var == "nx":
                x_label, y_label, z_label, zero_label = "mx", "my", "ny", "nx"
            elif eq_zero_var == "ny":
                x_label, y_label, z_label, zero_label = "mx", "my", "nx", "ny"
        elif eq_zero_var in ["px", "py", "sx", "sy"]:
            if eq_zero_var == "px":
                x_label, y_label, z_label, zero_label = "sx", "sy", "py", "px"
            elif eq_zero_var == "py":
                x_label, y_label, z_label, zero_label = "sx", "sy", "px", "py"
            elif eq_zero_var == "sx":
                x_label, y_label, z_label, zero_label = "px", "py", "sy", "sx"
            elif eq_zero_var == "sy":
                x_label, y_label, z_label, zero_label = "px", "py", "sx", "sy"
        return x_label, y_label, z_label, zero_label

    x1_label, y1_label, z1_label, zero1_label = get_labels(eq_zero_var1)
    x2_label, y2_label, z2_label, zero2_label = get_labels(eq_zero_var2)
    return dict(x1l=x1_label, y1l=y1_label, z1l=z1_label, zero1l=zero1_label, x2l=x2_label, y2l=y2_label, z2l=z2_label,
                zero2l=zero2_label)


@app.callback(
    Output('memory_table_data_and_styles', 'data'),
    Input('id_graph1', 'clickData'), Input('id_graph2', 'clickData'),
    State("id_table_show_sols0", 'style_cell'), State("id_table_show_sols0", 'style_header'),
    State("id_table_show_sols0", 'data'), State('memory_axis_labels', 'data'),
    Input('memory_input_panel_components_values', 'data'),
    prevent_initial_call=True)
def store_table_data_and_styles(cd1, cd2, table_cells_s, table_header_s, table_data, axis_labels, _input_parameters):
    printinfo("store_table_data_and_styles", ctx.triggered_id)

    # Store table data and table cells and header styles
    # When stored data and styles are updated then GUI table is updated (check: 'manage_content_style_in_tab1')

    # local function to update old table data to new by clicked point
    def update_table(old_data, cd, x_label, y_label, z_label, eq_zero_var):
        # The function is called:
        #  - every time a point as clicked,
        #  - when the graph is first time crated
        # When the graph is first time crated the 'click data' is None, so condition "if cd is None" is necessary to
        # executing without errors. If condition 'cd = None" is true, then table data is set to empty (as '-').
        # Otherwise table is updated.
        if cd is None:
            return [dict(**{name: '-' for name in table_header_names})]
        else:
            new_data = old_data
            x, y, z = cd["points"][0]["x"], cd["points"][0]["y"], cd["points"][0]["z"]
            new_data[0][x_label] = round(x, 1)
            new_data[0][y_label] = round(y, 1)
            new_data[0][z_label] = round(z, 1)
            new_data[0][eq_zero_var] = 0
            return new_data

    # if function is called by changed input parameters from input panel the table is set to no-active mode
    if ctx.triggered_id == 'memory_input_panel_components_values':
        return dict(cell=table_cells_s, header=table_header_s,
                    data=[dict(**{name: '-' for name in table_header_names})])
    # if function is called by 'click data' the table is updating to new values
    elif ctx.triggered_id == 'id_graph1':
        return dict(cell=table_cells_s, header=table_header_s,
                    data=update_table(table_data, cd1, axis_labels["x1l"], axis_labels["y1l"],
                                      axis_labels["z1l"], axis_labels["zero1l"]))
    elif ctx.triggered_id == 'id_graph2':
        return dict(cell=table_cells_s, header=table_header_s,
                    data=update_table(table_data, cd2, axis_labels["x2l"], axis_labels["y2l"],
                                      axis_labels["z2l"], axis_labels["zero2l"]))


@app.callback(
    Output('memory_input_panel_components_values', 'data'),
    Input('id_dd1', 'value'), Input('id_dd2', 'value'),
    Input('id_AB', 'value'), Input('id_BC', 'value'), Input('id_CD', 'value'),
    prevent_initial_call=True)
def store_input_panel_components_values(dd1, dd2, ab, bc, cd):
    printinfo("store_input_panel_components_values", ctx.triggered_id)

    # Store input parameters from input panel
    # Stored data is used to:
    #  - markers style change in 'id_graph1' and 'id_graph2' after changed the input parameters into
    #    input panel (check: 'select_active_points')
    #  - change the table data, cells and header styles to no-active after changed the input parameters into
    #    input panel (check: 'store_table_data_and_styles')
    return dict(ab_val=ab, bc_val=bc, cd_val=cd, dd1_val=dd1, dd2_val=dd2)


@app.callback(
    Output("memory_input_panel_components_styles", "data"),
    Input('cont-tabs', 'value'), Input('id_AB', 'disabled'), Input('id_BC', 'disabled'), Input('id_CD', 'disabled'),
    Input('id_dd1', 'disabled'), Input('id_dd2', 'disabled'), Input('id_dropdown1_div', 'className'),
    Input('id_dropdown2_div', 'className'), Input('id_btn', 'disabled'), Input("id_table_show_sols0", 'style_cell'),
    Input("id_table_show_sols0", 'style_header'), Input("id_table_show_sols0", 'data'), Input('cont-tab3', 'disabled'))
def store_input_panel_components_styles(tab_name, ab_s, bc_s, cd_s, dd1_s, dd2_s, dd1_cn, dd2_cn, calc_btn_s,
                                        table_cells_s, table_header_s, tab_data, tab_model_params_disabled):
    printinfo("store_input_panel_components_styles", ctx.triggered_id)

    # Store components styles from input panel when sub-tab is switched
    # The task of the function is to save a part of view. The view will load when tab "choose solutions"
    # are switched again (check: 'manage_content_style_in_tab1')
    if tab_name == 'choose_solutions':
        return dict(ab_style=ab_s, bc_style=bc_s, cd_style=cd_s, dd1_style=dd1_s, dd2_style=dd2_s, dd1_classname=dd1_cn,
                    dd2_classname=dd2_cn, calc_btn_style=calc_btn_s, conttab3dis=tab_model_params_disabled,
                    tabcell_style=table_cells_s, tabheader_style=table_header_s, tab_data=tab_data)
    else:
        return no_update


if __name__ == "__main__":
    app.run_server(port=4050, debug=False)
    # app.run_server(host='127.0.0.1', port='7080', proxy=None, debug=True, dev_tools_ui=True,
    #                dev_tools_props_check=None, dev_tools_serve_dev_bundles=None, dev_tools_hot_reload=None,
    #                dev_tools_hot_reload_interval=None, dev_tools_hot_reload_watch_interval=None,
    #                dev_tools_hot_reload_max_retry=None, dev_tools_silence_routes_logging=None,
    #                dev_tools_prune_errors=None)
