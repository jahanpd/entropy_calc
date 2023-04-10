use sqlx::{
    sqlite::{SqliteConnectOptions, SqliteJournalMode, SqlitePoolOptions, SqliteSynchronous},
    Pool, Sqlite,
};
use core::time::Duration;
use std::{str::FromStr, hash::Hash};
use std::sync::mpsc::{Receiver, Sender};
use futures::executor::block_on;
use bincode::{serialize, deserialize};

#[derive(Debug)]
pub struct EntropyData {
    pub ent: f64,
    pub celltype: String,
    pub gene: String,
    pub time: String,
    // pub sequence: Vec<Vec<i64>>
}

// output database object
pub struct EntropyDB {
    path: String,
    url: String,
    connections: u32,
    timeout: Duration,
    pool: Pool<Sqlite>,
    receiver: Receiver<Option<EntropyData>>
}

impl EntropyDB {
    pub async fn new(
        path: &str,
        connections: u32,
        timeout: u64,
        receiver: Receiver<Option<EntropyData>>
    ) -> Result<EntropyDB, Box<dyn std::error::Error>> {
        let url = format!("sqlite://{}", path);

    let connection_options = SqliteConnectOptions::from_str(&url)?
        .create_if_missing(true)
        .journal_mode(SqliteJournalMode::Wal)
        .synchronous(SqliteSynchronous::Off)
        .busy_timeout(Duration::from_secs(timeout));

    let sqlite_pool = SqlitePoolOptions::new()
        .max_connections(connections)
        .idle_timeout(Duration::from_secs(timeout))
        .connect_with(connection_options)
        .await?;

    sqlx::query("CREATE TABLE IF NOT EXISTS entropy (id TEXT, entropy REAL, gene TEXT, cell TEXT, time TEXT)")
        .execute(&sqlite_pool)
        .await?;
    sqlx::query("pragma temp_store = memory;")
        .execute(&sqlite_pool)
        .await?;
    sqlx::query("pragma mmap_size = 30000000000;")
        .execute(&sqlite_pool)
        .await?;
    sqlx::query("pragma page_size = 4096;")
        .execute(&sqlite_pool)
        .await?;

    // TODO create table migration

    return Ok(EntropyDB {
        path: String::from(path),
        url: String::from(url),
        connections,
        timeout: Duration::from_secs(timeout),
        pool: sqlite_pool,
        receiver
        })
    }

    pub async fn insert(&self, data: i32) {

    }

    pub fn listen(&self) {
        let rt = tokio::runtime::Builder::new_multi_thread()
            .worker_threads(12)
            .enable_all()
            .build().unwrap();
        while let Ok(data) = self.receiver.recv() {
            if data.is_none() {break};
            // dbg!(data);
            let d = data.unwrap();
            let id = vec![d.celltype.clone(), d.gene.clone()].join(":");
            // let seq: Vec<u8> = serialize(&d.sequence).unwrap();
            let _ = rt.block_on(sqlx::query(
                "INSERT OR REPLACE INTO entropy (id, entropy, gene, cell, time) VALUES (?, ?, ?, ?, ?)",
            )
            .bind(id)
            .bind(d.ent)
            .bind(d.gene)
            .bind(d.celltype)
            .bind(d.time)
            .execute(&self.pool)
            ).unwrap();
        }
    }
}
