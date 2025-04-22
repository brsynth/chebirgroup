from sqlalchemy import and_, Boolean, Column, create_engine, Float, Integer, String
from sqlalchemy.ext.declarative import declarative_base


Base = declarative_base()


# Define the Compound model
class Compound(Base):
    __tablename__ = "compound"

    cid = Column(Integer, primary_key=True, autoincrement=True)
    inchi = Column(String, nullable=True)
    is_biochemical = Column(Boolean, nullable=True)
    exact_mol_wt = Column(Float, nullable=True)

    def __repr__(self):
        return f"<Compound(id='{self.id}', cid='{self.cid}', inchi='{self.inchi}', is_biochemical='{self.is_biochemical}', exact_mol_wt={self.exact_mol_wt})>"
