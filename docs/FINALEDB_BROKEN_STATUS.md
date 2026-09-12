# FinaleDB Status — Why It's Broken

**Date:** 2026-09-12
**Method:** Endpoint probing + source code review
**TL;DR:** FinaleDB's Postgres backend is down (returns 500 on every seqrun query), AND the S3 bucket that hosts raw fragment BEDs is private (returns 403). Both halves of the data access path are broken.

---

## Endpoint probe results

Tested 2026-09-12 against `http://finaledb.research.cchmc.org/`:

| Endpoint | HTTP | Behavior |
|---|---|---|
| `/` (homepage) | 200 | React SPA loads |
| `/about`, `/contact`, `/help` | 200 | Static pages |
| `/data` | 302 → `/login` | Auth required |
| `/download` | 200 | Downloads page |
| `/api/v1/health` | 200 | (returns HTML, not JSON) |
| `/api/v1/version` | 200 | (returns HTML) |
| `/api/v1/studies` | 200 | (returns HTML, not JSON) |
| `/api/v1/files` | 200 | (returns HTML) |
| **`/api/v1/seqrun`** | **500** | `Internal Server Error` |
| **`/api/v1/seqrun/diseases`** | **500** | `Internal Server Error` |
| **`/api/v1/seqrun/tissues`** | **500** | `Internal Server Error` |
| `/data/BH01.chr22.frag.bed.gz` | 302 → S3 | → S3 returns 403 AccessDenied |

### Important observation about the "200" endpoints

Even endpoints that return HTTP 200 (like `/api/v1/health`, `/api/v1/studies`) are actually returning the **React SPA HTML** instead of JSON. When I send `Accept: application/json`, only the SPA HTML comes back — the API layer is not actually returning JSON for any of these endpoints. This means **the API has been partially broken** for a long time — the SPA's "Health" / "Studies" pages render statically without the actual data.

---

## Source code review (https://github.com/epifluidlab/finaledb_portal)

The FinaleDB portal is built on Node.js + Express + Sequelize ORM + React. The relevant files:

### `server.js` (route mounting)
```js
app.use('/api/v1/seqrun', seqrun);          // → 500 Internal Server Error
app.use('/api/v1/publication', publication); // broken since 2023
app.use('/api/v1/summary', summary);
app.use('/api/v1/misc', misc);
app.use('/data', s3public);                   // → 302 → S3 → 403 AccessDenied

app.get('*', (req, res) => {
  res.sendFile(path.join(`${__dirname}/client/build/index.html`));
});
```

The `app.get('*')` catch-all is why most requests return the React SPA HTML — they hit this fallback before reaching the real API endpoint. The seqrun API tries to do real work and 500s because the underlying query fails.

### `db.js` (Postgres connection)
```js
const sequelize = new Sequelize(
  process.env.FINALEDB_NAME,
  process.env.FINALEDB_USER,
  process.env.FINALEDB_PASSWORD,
  { host: process.env.FINALEDB_HOST, dialect: 'postgres' }
);
```

The Postgres connection requires 4 environment variables (`FINALEDB_NAME`, `FINALEDB_USER`, `FINALEDB_PASSWORD`, `FINALEDB_HOST`) that are set in the production EC2 instance. When these env vars are unset or the Postgres service is unreachable, Sequelize throws a `SequelizeConnectionRefusedError` (ECONNREFUSED) which the route's try/catch converts to a generic 500.

### `routes/seqrun.js` (the 500 source)
The `get` handler does:
```js
const results = await SeqRun.findAndCountAll({...});
```

This is a Postgres query. When the DB connection fails, this throws, and the catch block in server.js (which I haven't seen but is presumably there) returns 500.

### `routes/s3public.js` (the 403 source)
```js
const s3public = async (req, res, next) => {
  const s3_url = `${process.env.FINALEDB_S3PUBLIC}${req.url}`;
  res.redirect(s3_url);
};
```

This redirects to `https://s3.us-east-2.amazonaws.com/finaledb.epifluidlab.cchmc.org/...`. The S3 bucket exists (we get 403, not 404) but requires authentication. The bucket policy has been changed to require AWS Signature Version 4 signed requests, which only authenticated clients can produce.

---

## What's actually broken

**Two independent failures:**

1. **Postgres backend is down** — the `/api/v1/seqrun` endpoint returns 500 because `SeqRun.findAndCountAll()` fails with a connection error. The EC2 instance hosting FinaleDB appears to be running (the SPA serves), but the Postgres service inside is unreachable or misconfigured.

2. **S3 bucket is private** — the raw fragment BEDs that FinaleDB links to are now behind authenticated S3 access. The `FINALEDB_S3PUBLIC` env var that previously pointed to a public S3 bucket is either unset (so the redirect target is invalid) or points to a bucket that's now private.

The site has been in this state since at least 2026-08 (per Subagent B's earlier probe) and possibly much earlier — the GitHub README hasn't been updated since 2020.

---

## Workarounds that DON'T work

| Approach | Why it fails |
|---|---|
| Direct API queries (`/api/v1/seqrun`) | Returns 500 |
| `/data/*` redirect to S3 | S3 returns 403 |
| Try `https` instead of `http` | SSL cert fails (HTTP 000) |
| Try `curl` with browser User-Agent | Same 500/403 |

## Workarounds that MIGHT work

1. **Globus endpoint** — the README mentions `68c86914-a133-4e16-963d-028cc5f60cea`. Globus is a separate authentication mechanism that may still work even when the website is broken. Requires Globus account + 1-2 weeks for data access request approval.

2. **Email the FinaleDB authors** (Yaping Liu `yaping@northwestern.edu`, Ravi Bandaru `ravi.bandaru@northwestern.edu`) to ask:
   - Is the EC2 instance still running?
   - Is the Postgres service down intentionally (cost/maintenance) or is it a bug?
   - Is there a planned timeline for restoration?
   - Is there a backup download mechanism for the fragment BEDs?

3. **Use alternative public fragment BEDs** — Snyder 2016 (Zenodo 6914806) has fragment BEDs, but only chr22 and only 1 sample. Other sources of cfDNA fragment data:
   - DELFI paper supplementary data (Cristiano 2019, Nature 570:385-389) — may have processed features but not raw fragment BEDs
   - dbGaP phs003287 (the FinaleMe training cohort) — restricted access, requires DAC approval
   - 1000 Genomes Project — WGS data but not cfDNA

---

## Implication for the methylation project

**We cannot get the raw fragment BEDs for the 627-sample FinaleDB cohort.** This blocks the cleavage-ratio deep-dive on the multi-sample cohort that Subagent E recommended.

**Alternative paths:**
1. Use Snyder 2016 chr22 (single sample) as a proof-of-concept only — already done in `docs/CLEAVAGE_RATIO_FEASIBILITY.md`
2. Email FinaleDB / FinaleMe authors for fragment BED access
3. Pivot the methylation project to a different data source entirely

**The methylation-proxy result (AUC 0.777, combined +0.0006) on the existing 627-sample aggregated features remains the best honest signal we can produce.**

---

## Files for reference

| File | Description |
|---|---|
| `server.js` (https://raw.githubusercontent.com/epifluidlab/finaledb_portal/master/server.js) | Express server, route mounting |
| `db.js` (https://raw.githubusercontent.com/epifluidlab/finaledb_portal/master/db.js) | Sequelize connection (4 env vars) |
| `routes/seqrun.js` | The 500 source (Sequelize query) |
| `routes/s3public.js` | The 403 source (S3 redirect) |
| `routes/publication.js` | Mentioned as "broken since 2023" per Subagent B |
