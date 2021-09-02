/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#include"DB.h"
#include<cstring>
using namespace DB;
HParams ParamsDB::fetch(const char* stock_name){
    sqlite3_stmt *res_handle;
    HParams ret_val;
    int rc;
    if((rc=sqlite3_prepare_v2(db, fetch_row, -1, &res_handle, 0))==SQLITE_OK) sqlite3_bind_text(res_handle, sqlite3_bind_parameter_index(res_handle, "@name"), stock_name,strlen(stock_name),nullptr);
    else throw DBError(db,sqlite3_errmsg(db),rc);
    if((rc=sqlite3_step(res_handle))==SQLITE_ROW)
        ret_val={sqlite3_column_double(res_handle,V0),
                 sqlite3_column_double(res_handle,VM),
                 sqlite3_column_double(res_handle,RHO),
                 sqlite3_column_double(res_handle,KAPPA),
                 sqlite3_column_double(res_handle,SIGMA)};
    else if(rc==SQLITE_DONE) throw DBException();
    else throw DBError(db,sqlite3_errmsg(db),rc);
    if((rc=sqlite3_finalize(res_handle))!=SQLITE_OK) throw DBError(db,sqlite3_errmsg(db),rc);
    return ret_val;
}
void ParamsDB::insertupdate(HParams* p, const char* stock_name){
    sqlite3_stmt *res_handle;
    int rc;
    if((rc=sqlite3_prepare_v2(db, insert_update, -1, &res_handle, 0))==SQLITE_OK){
        sqlite3_bind_text(res_handle,   sqlite3_bind_parameter_index(res_handle, "@name"),stock_name,strlen(stock_name),nullptr);
        sqlite3_bind_double(res_handle, sqlite3_bind_parameter_index(res_handle, "@v0"), p->v_0);
        sqlite3_bind_double(res_handle, sqlite3_bind_parameter_index(res_handle, "@vm"), p->v_m);
        sqlite3_bind_double(res_handle, sqlite3_bind_parameter_index(res_handle, "@rho"), p->rho);
        sqlite3_bind_double(res_handle, sqlite3_bind_parameter_index(res_handle, "@kappa"), p->kappa);
        sqlite3_bind_double(res_handle, sqlite3_bind_parameter_index(res_handle, "@sigma"), p->sigma);
    }else throw DBError(db,sqlite3_errmsg(db),rc);
    if((rc=sqlite3_step(res_handle))!=SQLITE_DONE||(rc=sqlite3_finalize(res_handle))!=SQLITE_OK) throw DBError(db,sqlite3_errmsg(db),rc);
}
/*
void ParamsDB::update(HParams* p, const char* stock_name){
    insertupdate(p,stock_name,update_row);
}
void ParamsDB::insert(HParams* p, const char* stock_name){
    insertupdate(p,stock_name,insert_row);
}*/
