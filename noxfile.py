import nox

@nox.session(python=["3.9", "3.10"])
def tests(session):
    session.run("poetry", "install", external=True)
    session.run("pytest")
