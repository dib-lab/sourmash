import pytest


@pytest.fixture(params=[True, False])
def cli(request):
    return request.param


@pytest.fixture(params=[True, False])
def track_abundance(request):
    return request.param


@pytest.fixture(params=[True, False])
def dayhoff(request):
    return request.param


@pytest.fixture(params=[True, False])
def hp(request):
    return request.param


@pytest.fixture(params=[2, 5, 10])
def n_children(request):
    return request.param


# --- BEGIN - Only run tests using a particular fixture --- #
# Cribbed from: http://pythontesting.net/framework/pytest/pytest-run-tests-using-particular-fixture/
def pytest_collection_modifyitems(items, config):
    fixture_name = config.option.usesfixture
    if fixture_name is not None:
        selected_items = []
        deselected_items = []

        for item in items:
            if fixture_name in getattr(item, 'fixturenames', ()):
                selected_items.append(item)
            else:
                deselected_items.append(item)
        config.hook.pytest_deselected(items=deselected_items)
        items[:] = selected_items

def pytest_addoption(parser):
    parser.addoption("--usesfixture",
                     action="store",
                     default=None,
                     help="just run tests that use a particular fixture")
# --- END - Only run tests using a particular fixture --- #
