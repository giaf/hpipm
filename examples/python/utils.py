
import sys
import os


def check_env():
    """ Check that env.sh has been run."""

    env_run = os.getenv('ENV_RUN')
    if env_run != 'true':
        print('ERROR: env.sh has not been sourced! Before executing this example, run:')
        print('source env.sh')
        sys.exit(1)
    return


def is_travis_run() -> bool:
    travis_env_var = os.getenv('TRAVIS_RUN')
    return travis_env_var == 'true'