-- MR GUI project database, schema v1 (contracts C8, C9)
-- Per-connection pragmas set by the engine: journal_mode=WAL, synchronous=FULL, foreign_keys=ON.

CREATE TABLE meta (
  key   TEXT PRIMARY KEY,
  value TEXT NOT NULL
);
INSERT INTO meta(key, value) VALUES ('schema_version', '1'), ('contract_version', '1.0.0');

CREATE TABLE blobs (
  sha256     TEXT PRIMARY KEY CHECK (length(sha256) = 64),
  size       INTEGER NOT NULL CHECK (size >= 0),
  media_type TEXT NOT NULL,
  created_at TEXT NOT NULL
);

CREATE TABLE datasets (
  dataset_id         TEXT PRIMARY KEY,
  card_json          TEXT NOT NULL,
  variant_table_blob TEXT REFERENCES blobs(sha256),
  created_at         TEXT NOT NULL
);

CREATE TABLE instrument_sets (
  instrument_set_id TEXT PRIMARY KEY,
  card_json         TEXT NOT NULL,
  variants_blob     TEXT REFERENCES blobs(sha256),
  created_at        TEXT NOT NULL
);

CREATE TABLE tasks (
  task_id         TEXT PRIMARY KEY,
  task_type       TEXT NOT NULL CHECK (task_type IN ('download','import','ancestry','instruments','harmonise',
                    'ld_matrix','ldsc','check','method','figure','table','report')),
  state           TEXT NOT NULL CHECK (state IN ('pending','running','done','failed','blocked','cancelled')),
  fingerprint     TEXT CHECK (fingerprint IS NULL OR length(fingerprint) = 64),
  method_id       TEXT,
  adapter_version TEXT,
  seed            INTEGER,
  attempt         INTEGER NOT NULL DEFAULT 0,
  started_at      TEXT,
  finished_at     TEXT,
  updated_at      TEXT NOT NULL
);
CREATE INDEX tasks_state ON tasks(state);

CREATE TABLE task_deps (
  task_id    TEXT NOT NULL REFERENCES tasks(task_id) ON DELETE CASCADE,
  depends_on TEXT NOT NULL REFERENCES tasks(task_id) ON DELETE CASCADE,
  PRIMARY KEY (task_id, depends_on),
  CHECK (task_id <> depends_on)
);

CREATE TABLE results (
  task_id            TEXT PRIMARY KEY REFERENCES tasks(task_id) ON DELETE CASCADE,
  fingerprint        TEXT NOT NULL CHECK (length(fingerprint) = 64),
  status             TEXT NOT NULL CHECK (status IN ('ok','warn','failed','blocked')),
  result_blob        TEXT NOT NULL REFERENCES blobs(sha256),
  variant_level_blob TEXT REFERENCES blobs(sha256),
  created_at         TEXT NOT NULL
);

CREATE TABLE ui_state (
  key        TEXT PRIMARY KEY,
  value_json TEXT NOT NULL,
  updated_at TEXT NOT NULL
);

CREATE TABLE audit_log (
  id          INTEGER PRIMARY KEY AUTOINCREMENT,
  at          TEXT NOT NULL,
  event       TEXT NOT NULL,
  detail_json TEXT NOT NULL DEFAULT '{}'
);
