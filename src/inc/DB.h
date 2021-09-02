/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
#pragma once
#include <sqlite3.h>
#include"HDistribution.h"
namespace DB{
    struct DBException : public std::exception
    {
        const char * what () const throw ()
        {
            //TODO copy error_msg and sqlite3_free error_msg
            return "try again";
        }
    };
    struct DBError : public std::exception
    {
        const char* m_error_msg;
        DBError(sqlite3 *db,const char* error_msg,int rc){
            sqlite3_close(db);
            m_error_msg=error_msg;
            std::cout<<"Error: "<<rc<<std::endl;
        }
        const char * what () const throw ()
        {
            //TODO copy error_msg and sqlite3_free error_msg
            return m_error_msg;
        }
    };
    //const std::vector<std::string> names={"Name","V0","Vm","Rho","Kappa","Sigma","DaysToExpiry"};
    const char create_table[]="CREATE TABLE IF NOT EXISTS HestonParameters (Name TEXT NOT NULL UNIQUE, V0 REAL NOT NULL, Vm REAL NOT NULL, Rho REAL NOT NULL, Kappa REAL NOT NULL, Sigma REAL NOT NULL, ComputedAt INTEGER NOT NULL);";
    const char fetch_row[]="SELECT * FROM HestonParameters WHERE Name=@name";
    const char insert_update[]="INSERT INTO HestonParameters VALUES(@name,@v0,@vm,@rho,@kappa,@sigma,strftime('%s','now')) ON CONFLICT(Name) DO UPDATE SET V0=@v0,Vm=@vm,Rho=@rho,Kappa=@kappa,Sigma=@sigma,ComputedAt=strftime('%s','now'); ";
    enum request{
        NAME=0,
        V0=1,
        VM=2,
        RHO=3,
        KAPPA=4,
        SIGMA=5,
        COMPUTEDAT=6
    };
    class ParamsDB{
        sqlite3 *db;
    public:
        ParamsDB(){
            char* err_msg=0;
            int rc;
            if((rc=sqlite3_open("HParams.db", &db))!=SQLITE_OK
                ||(rc=sqlite3_exec(db, create_table, 0, 0, &err_msg))!=SQLITE_OK) 
                throw DBError(db,err_msg,rc);
        };
        HParams fetch(const char* stock_name);
        void insertupdate(HParams* p, const char* stock_name);
        ~ParamsDB(){
            sqlite3_close(db);
        };
    };
}
