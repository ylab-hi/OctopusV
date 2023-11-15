import nox

@nox.session(python="3.9")
def tests_39(session):
    session.run("poetry", "install", external=True)
    session.run("pytest", "--cov", "--cov-report=xml")

@nox.session(python="3.10")
def tests_310(session):
    session.run("poetry", "install", external=True)
    session.run("pytest", "--cov", "--cov-report=xml")

@nox.session
def coverage(session):
    session.install("coverage")
    session.run("coverage", "html")