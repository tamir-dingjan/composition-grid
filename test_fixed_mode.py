from doctest import OutputChecker
import os
import shutil
import pandas as pd
import glob
import pytest
import composition_grid


def execute(config_file):
    data_path = os.path.join(os.getcwd(), "tests", "data.csv")
    config_path = os.path.join(os.getcwd(), "tests", config_file)
    output_path = os.path.join(os.getcwd(), "tests", os.path.basename(os.path.splitext(config_path)[0]))
    
    try:
        os.mkdir(output_path)
    except FileExistsError:
        shutil.rmtree(output_path)
        os.mkdir(output_path)

    composition_grid.generate_grids(output_path, data_path, config_path)
    return [x for x in glob.glob(os.path.join(output_path, "*csv"))]


def load_grid_composition_file(grid_filepath):
    # Test if output loadable
    try:
        grid = pd.read_csv(grid_filepath)
    except:
        raise
    return grid


def test_leaflet_sum():
    output_files = execute("config1.csv")

    for grid_filepath in output_files:
        # Run general tests on each file
        # Test if output generated
        assert os.path.isfile(grid_filepath)

        # Test if output loadable
        grid = load_grid_composition_file(grid_filepath)
        
        assert grid["inside"].sum() == pytest.approx(100.0)
        assert grid["outside"].sum() == pytest.approx(100.0)

        
@pytest.mark.parametrize(('file_index', 'conc'), [(0, 5.0),(1, 10.0),(2, 15.0),(3, 20.0),(4, 25.0)])
def test_inside_leaflet(file_index, conc):
    output_files = execute("config1.csv")
    grid = load_grid_composition_file(sorted(output_files)[file_index])
    assert grid["inside"].tolist()[0] == pytest.approx(conc)


@pytest.mark.parametrize(('file_index', 'conc'), [(0, 2.5),(1, 5.0),(2, 7.5),(3, 10.0),(4, 12.5)])
def test_outside_leaflet(file_index, conc):
    output_files = execute("config1.csv")
    grid = load_grid_composition_file(sorted(output_files)[file_index])
    assert grid["outside"].tolist()[0] == pytest.approx(conc)
