
#[derive(Debug, Clone)]
pub struct Gene {
    pub gene_id: String,
    pub contig: String,
    pub start: u64,
    pub end: u64
}

#[derive(Debug, Clone)]
pub struct StoreGenes {
    pub gene_id: Vec<String>,
    pub contig: Vec<String>,
    pub start: Vec<u64>,
    pub end: Vec<u64>
}

impl StoreGenes {
    pub fn new(csv: String) -> StoreGenes {
        let mut reader = csv::Reader::from_reader(csv.as_bytes());
        let mut genes = vec![];
        let mut contig = vec![];
        let mut start = vec![];
        let mut end = vec![];
        for record in reader.records() {
            let record = record.unwrap();
            genes.push(record[0].into());
            contig.push(record[1].into());
            start.push(record[2].parse::<u64>().unwrap());
            end.push(record[3].parse::<u64>().unwrap());
        }
        return StoreGenes {
            gene_id: genes,
            contig,
            start,
            end
        }
    }

    pub fn get_gene(&self, gene_id: String) -> Option<Gene> {
        match self.gene_id.iter().position(|r| r  == &gene_id) {
            Some(val) => return Some( Gene{
                gene_id: self.gene_id[val].clone(),
                contig: self.contig[val].clone(),
                start: self.start[val],
                end: self.end[val]
            } ),
            None => return None
        }
    }

    pub fn genes(&self) -> impl Iterator<Item = Gene> + '_ {
        return (0..self.gene_id.len()).map(|i| Gene {
            gene_id: self.gene_id[i].clone(),
            contig: self.contig[i].clone(),
            start: self.start[i],
            end: self.end[i],
        });
    }
}
