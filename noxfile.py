import nox

@nox.session(python="3.9")
def tests_39(session):
    session.run("poetry", "install", external=True)
    session.run("pytest", "--cov", "--cov-report=html")

@nox.session(python="3.10")
def tests_310(session):
    session.run("poetry", "install", external=True)
    session.run("pytest", "--cov", "--cov-report=html")
